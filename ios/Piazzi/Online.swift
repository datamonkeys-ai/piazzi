import Foundation

/// Public services used by the app (the same ones astromcp queries): IMCCE SkyBoT and VizieR.
enum Online {
    struct KnownAsteroid: Identifiable { var id: String { name }; let number: String?, name: String, ra: Double, dec: Double, vmag: Double, sep: Double, cls: String }

    static let vizierHosts = ["vizier.idia.ac.za", "vizier.cds.unistra.fr", "vizier.cfa.harvard.edu"]

    /// GET with retries (IRSA sometimes returns 503); rejects HTML/short responses when `minBytes` is given.
    static func get(_ url: URL, minBytes: Int = 0, tries: Int = 3, timeout: TimeInterval = 60) async throws -> Data {
        var last: Error = URLError(.badServerResponse)
        for t in 0..<tries {
            do {
                try Task.checkCancellation()
                var req = URLRequest(url: url); req.timeoutInterval = timeout
                let (d, r) = try await URLSession.shared.data(for: req)
                if let h = r as? HTTPURLResponse, h.statusCode == 200, d.count >= minBytes { return d }
                last = URLError(.badServerResponse)
            } catch { if error is CancellationError || (error as? URLError)?.code == .cancelled { throw error }; last = error }
            if t + 1 < tries { try await Task.sleep(nanoseconds: UInt64(2_000_000_000 * (t + 1))) }   // IRSA 503: increasing backoff
        }
        throw last
    }

    /// Builds an HTTPS URL with encoded parameters (URLComponents/URLQueryItem); "+" becomes %2B so it is not read
    /// as a space by the server. Throws URLError(.badURL) instead of crashing the app with a force-unwrap.
    static func makeURL(host: String, path: String, _ items: [(String, String)]) throws -> URL {
        var c = URLComponents()
        c.scheme = "https"; c.host = host; c.path = path
        c.queryItems = items.map { URLQueryItem(name: $0.0, value: $0.1) }
        c.percentEncodedQuery = c.percentEncodedQuery?.replacingOccurrences(of: "+", with: "%2B")
        guard let u = c.url else { throw URLError(.badURL) }
        return u
    }

    // ---------------------------------------------------------------- SkyBoT
    static func skybot(jd: Double, ra: Double, dec: Double, radiusDeg: Double, observer: String) async throws -> [KnownAsteroid] {
        let obs = observer.isEmpty || observer == "XXX" ? "500" : observer
        // observatory code comes from the user: sent encoded as a parameter, never interpolated into the URL
        let u = try makeURL(host: "ssp.imcce.fr", path: "/webservices/skybot/api/conesearch.php",
                            [("-ep", "\(jd)"), ("-ra", "\(ra)"), ("-dec", "\(dec)"), ("-rd", "\(radiusDeg)"), ("-mime", "text"),
                             ("-output", "object"), ("-observer", obs), ("-objFilter", "110"), ("-from", "Piazzi")])
        let text = String(decoding: try await get(u, tries: 2, timeout: 25), as: UTF8.self)
        var out: [KnownAsteroid] = []
        for l in text.split(separator: "\n") where !l.hasPrefix("#") && l.contains("|") {
            let c = l.split(separator: "|", omittingEmptySubsequences: false).map { $0.trimmingCharacters(in: .whitespaces) }
            guard c.count >= 8, let v = Double(c[5]), let sep = Double(c[7]) else { continue }
            out.append(KnownAsteroid(number: c[0] == "-" || c[0].isEmpty ? nil : c[0], name: c[1], ra: sexa(c[2]) * 15, dec: sexa(c[3]), vmag: v, sep: sep, cls: c[4]))
        }
        return out.sorted { $0.sep < $1.sep }
    }

    static func sexa(_ s: String) -> Double {
        let p = s.split(separator: " ").compactMap { Double($0) }
        guard p.count == 3 else { return Double(s) ?? 0 }
        let sg: Double = s.trimmingCharacters(in: .whitespaces).hasPrefix("-") ? -1 : 1
        return sg * (abs(p[0]) + p[1] / 60 + p[2] / 3600)
    }
}

// MARK: - "What is this?" (like astromcp's whats_here): SIMBAD + VSX + Gaia DR3 + SkyBoT in parallel
extension Online {
    struct Inspection {
        var simbad: [(id: String, type: String, sep: Double)] = []
        var vsx: [(name: String, type: String, range: String, period: String)] = []
        var gaia: (g: Double, bprp: Double?, plx: Double?, sep: Double)?
        var asteroid: KnownAsteroid?
        var failed: [String] = []
        var verdict: (title: String, detail: String, kind: Kind) {
            if let a = asteroid { return ((a.number.map { "(\($0)) " } ?? "") + a.name, L("Asteroid · %@ · V %@", a.cls, String(format: "%.1f", a.vmag)), .asteroid) }
            if let v = vsx.first { return (v.name, L("Variable star") + " \(v.type) · \(v.range)\(v.period.isEmpty ? "" : " · P = \(v.period) d")", .variable) }
            if let s = simbad.first { return (s.id, SimbadType.describe(s.type), s.type.contains("*") || s.type.lowercased() == "star" ? .star : .other) }
            if let g = gaia { return (L("Gaia DR3 star"), String(format: "G %.2f", g.g) + (g.bprp.map { String(format: " · BP−RP %.2f", $0) } ?? "") + (g.plx.flatMap { $0 > 0 ? L(" · ~%.0f light-years", 3261.6 / $0) : nil } ?? ""), .star) }
            return (L("Nothing known here"), L("No catalogued object within 6″ — noise, artifact or something new."), .unknown)
        }
        enum Kind { case asteroid, variable, star, other, unknown }
    }

    static let simbadHosts = ["simbad.cfa.harvard.edu", "simbad.cds.unistra.fr"]

    static func whatIsHere(ra: Double, dec: Double, jd: Double, observer: String, radiusArcsec: Double = 6) async -> Inspection {
        var r = Inspection()
        async let sb: [(String, String, Double)]? = simbad(ra: ra, dec: dec, radius: radiusArcsec)
        async let vx: [(String, String, String, String)]? = vsx(ra: ra, dec: dec, radius: max(radiusArcsec, 8))
        async let ga: (Double, Double?, Double?, Double)?? = gaia(ra: ra, dec: dec, radius: radiusArcsec)
        async let ast: [KnownAsteroid]? = try? skybot(jd: jd, ra: ra, dec: dec, radiusDeg: max(radiusArcsec, 10) / 3600, observer: observer)
        if let v = await sb { r.simbad = v.map { (id: $0.0, type: $0.1, sep: $0.2) } } else { r.failed.append("SIMBAD") }
        if let v = await vx { r.vsx = v.map { (name: $0.0, type: $0.1, range: $0.2, period: $0.3) } } else { r.failed.append("VSX") }
        if let v = await ga { r.gaia = v.map { (g: $0.0, bprp: $0.1, plx: $0.2, sep: $0.3) } } else { r.failed.append("Gaia") }
        if let v = await ast { r.asteroid = v.first } else { r.failed.append("SkyBoT") }
        return r
    }

    static func simbad(ra: Double, dec: Double, radius: Double) async -> [(String, String, Double)]? {
        let q = "SELECT TOP 3 main_id, otype, DISTANCE(POINT('ICRS',ra,dec),POINT('ICRS',\(ra),\(dec)))*3600 AS d FROM basic WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',\(ra),\(dec),\(radius / 3600)))=1 ORDER BY d"
        for h in simbadHosts {
            guard let u = try? makeURL(host: h, path: "/simbad/sim-tap/sync", [("REQUEST", "doQuery"), ("LANG", "ADQL"), ("FORMAT", "csv"), ("QUERY", q)]),
                  let d = try? await get(u, tries: 1, timeout: 15) else { continue }
            let lines = String(decoding: d, as: UTF8.self).split(separator: "\n").dropFirst()
            return lines.compactMap { l in
                let c = l.split(separator: ",", omittingEmptySubsequences: false).map { $0.trimmingCharacters(in: CharacterSet(charactersIn: "\" ")) }
                guard c.count >= 3, let s = Double(c[2]) else { return nil }
                return (c[0].replacingOccurrences(of: "  ", with: " "), c[1], s)
            }
        }
        return nil
    }

    static func vizierRows(_ source: String, ra: Double, dec: Double, radius: Double, cols: String) async -> [[String]]? {
        for h in vizierHosts {
            let c = "\(ra) \(dec >= 0 ? "+" : "")\(dec)"
            guard let u = try? makeURL(host: h, path: "/viz-bin/asu-tsv", [("-source", source), ("-c", c), ("-c.rs", "\(radius)"), ("-out", cols),
                                                                        ("-oc.form", "d"), ("-sort", "_r")]),
                  let d = try? await get(u, tries: 1, timeout: 15) else { continue }
            let lines = String(decoding: d, as: UTF8.self).split(separator: "\n", omittingEmptySubsequences: false).filter { !$0.hasPrefix("#") && !$0.isEmpty }
            guard lines.count >= 1 else { continue }
            return lines.dropFirst(3).map { $0.split(separator: "\t", omittingEmptySubsequences: false).map { $0.trimmingCharacters(in: .whitespaces) } }
        }
        return nil
    }

    static func vsx(ra: Double, dec: Double, radius: Double) async -> [(String, String, String, String)]? {
        guard let rows = await vizierRows("B/vsx/vsx", ra: ra, dec: dec, radius: radius, cols: "Name,Type,max,min,Period") else { return nil }
        return rows.filter { $0.count >= 5 && !$0[0].isEmpty }.map { r in
            let p = Double(r[4]).map { String(format: $0 < 10 ? "%.4f" : "%.1f", $0) } ?? ""
            return (r[0], r[1], "\(r[2])–\(r[3]) mag", p)
        }
    }

    static func gaia(ra: Double, dec: Double, radius: Double) async -> (Double, Double?, Double?, Double)?? {
        guard let rows = await vizierRows("I/355/gaiadr3", ra: ra, dec: dec, radius: radius, cols: "_r,Gmag,BP-RP,Plx") else { return nil }
        guard let r = rows.first(where: { $0.count >= 4 && Double($0[1]) != nil }) else { return .some(nil) }
        return .some((Double(r[1])!, Double(r[2]), Double(r[3]), Double(r[0]) ?? 0))
    }
}

/// Most common SIMBAD object types in Portuguese (localized labels).
enum SimbadType {
    static let map: [String: String] = ["*": "Star", "s*r": "Red supergiant", "G": "Galaxy", "QSO": "Quasar", "PN": "Planetary nebula",
        "HII": "HII region", "Cl*": "Star cluster", "OpC": "Open cluster", "GlC": "Globular cluster", "V*": "Variable star",
        "EB*": "Eclipsing binary", "**": "Double star", "WD*": "White dwarf", "SN": "Supernova", "Rad": "Radio source", "X": "X-ray source",
        "IR": "Infrared source", "PM*": "High proper-motion star", "RR*": "RR Lyrae variable", "Mi*": "Mira variable"]
    static func describe(_ t: String) -> String { "SIMBAD · " + (map[t].map { L($0) } ?? t) }
}

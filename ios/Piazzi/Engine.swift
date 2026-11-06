import Foundation
import UIKit

struct FrameResult: Decodable, Identifiable {
    var id: String { name + time }
    let name: String, time: String, jd: Double, ndet: Int, nref: Int, nmatch: Int, dra: Double?, dde: Double?, zp: Double?, nphot: Int, seeing: Double?, w: Int, h: Int // dra/dde/zp/seeing: null when not finite
    let dets: [[Double]], refs: [[Double]]
}
struct Measure: Decodable {
    let x: Double, y: Double, ra: String, dec: String, mag: Double, rad: Double, decd: Double, jd: Double, time: String
    // fields added in 0.2.0 (optional: σ and trail length are null when undefined)
    let trksub: String?, sig_ra: Double?, sig_dec: Double?, trail_len: Double?, remarks: String?
}
struct Mover: Decodable, Identifiable {
    var id: String { label + String(rate ?? 0) }
    var desig: String; let trksub: String?; let rate: Double?, pa: Double?, meas: [Measure]
    // display label: without a designation (unidentified object) uses the track's temporary designation (trkSub)
    var label: String { !desig.isEmpty ? desig : (trksub.flatMap { $0.isEmpty ? nil : $0 } ?? L("object")) }
}
struct ReduceResult: Decodable { let catalog: Int?; let frames: [FrameResult]?; var movers: [Mover]?; let error: String?; let code: String?; let detail: String? }
extension ReduceResult { var errorText: String? { error == nil ? nil : coreError(code: code, message: error, detail: detail) } }

enum StepState: Equatable { case pending, running, done(String), failed(String) }
struct Step: Identifiable { let id: Int; let title: String; let icon: String; var state: StepState = .pending }
enum Phase { case home, running, results }

/// Equipment profile (telescope + camera + observatory).
struct Profile: Codable, Equatable {
    var name: String, focal: Double, pixelW: Double, pixelH: Double, flip: Int, mpc: String
    var summary: String { String(format: "%.0f mm · %.1f×%.1f µm · %.2f″/px · %@", focal, pixelW * 1000, pixelH * 1000, pixelW / focal * 206264.8, mpc) }
    static let builtIn: [Profile] = []
}

/// Bridge to the C++ core (detection, plate, photometry and motion) + online services.
@MainActor final class Engine: ObservableObject {
    private let s = astro_new()!
    @Published var images: [URL] = []
    @Published var catalog: URL?
    @Published var cfgName = "default"
    /// built-in example session: data is not the user's — MPC submission blocked
    @Published var isExample = false
    @Published var ra = ""
    @Published var dec = ""
    @Published var object = ""
    @Published var result: ReduceResult?
    @Published var known: [Online.KnownAsteroid] = []
    @Published var identified: [String: String] = [:]      // object desig -> identification description
    @Published var status = L("Import the images (.st6/.fit) and, if you like, a .cfg and the catalog.") {
        didSet { let a = ProcessInfo.processInfo.arguments; if a.contains("-gaia") { print("STATUS:", status) } }
    }
    @Published var busy = false
    private var work: Task<Void, Never>?
    @Published var phase: Phase = .home
    @Published var steps: [Step] = Engine.freshSteps()
    @Published var needCenter = false
    @Published var profileName = UserDefaults.standard.string(forKey: "profile") ?? "Default"
    @Published var profiles: [Profile] = {
        let saved = (try? JSONDecoder().decode([Profile].self, from: UserDefaults.standard.data(forKey: "profiles") ?? Data())) ?? []
        return Profile.builtIn + saved
    }()
    func apply(profile p: Profile) {
        set("focal_length_mm", String(p.focal)); set("pixel_width_mm", String(p.pixelW)); set("pixel_height_mm", String(p.pixelH))
        set("flip_horizontal", String(p.flip)); set("mpc_code", p.mpc); profileName = p.name; UserDefaults.standard.set(p.name, forKey: "profile")
    }
    func saveProfile(named n: String) {
        let p = Profile(name: n, focal: Double(get("focal_length_mm")) ?? 0, pixelW: Double(get("pixel_width_mm")) ?? 0, pixelH: Double(get("pixel_height_mm")) ?? 0,
                        flip: Int(get("flip_horizontal")) ?? 0, mpc: get("mpc_code"))
        profiles.removeAll { $0.name == n }; profiles.append(p); persistProfiles(); profileName = n; UserDefaults.standard.set(n, forKey: "profile")
    }
    func deleteProfiles(at idx: IndexSet) { profiles.remove(atOffsets: idx); profiles.insert(contentsOf: Profile.builtIn.filter { b in !profiles.contains(b) }, at: 0); persistProfiles() }
    private func persistProfiles() {
        let user = profiles.filter { p in !Profile.builtIn.contains(p) }
        UserDefaults.standard.set(try? JSONEncoder().encode(user), forKey: "profiles")
    }
    static func freshSteps() -> [Step] {
        [Step(id: 0, title: L("Images"), icon: "photo.stack"), Step(id: 1, title: L("Gaia DR3 catalog"), icon: "star.circle"),
         Step(id: 2, title: L("Astrometric plate"), icon: "scope"), Step(id: 3, title: L("Moving objects"), icon: "arrow.up.right.circle"),
         Step(id: 4, title: L("Identification"), icon: "checkmark.seal")]
    }
    private func step(_ i: Int, _ st: StepState) { steps[i].state = st }
    var frames: [FrameResult] { result?.frames ?? [] }
    var movers: [Mover] { result?.movers ?? [] }

    init() { if let p = profiles.first(where: { $0.name == profileName }) { apply(profile: p) } }
    deinit { astro_free(s) }

    func get(_ k: String) -> String { String(cString: astro_get(s, k)) }
    func set(_ k: String, _ v: String) { astro_set(s, k, v); objectWillChange.send() }

    func add(_ url: URL) {
        let ext = url.pathExtension.lowercased()
        if ext == "cfg" { astro_load_cfg(s, url.path); cfgName = url.lastPathComponent; objectWillChange.send() }
        else if ext == "txt" || ext == "tsv" { catalog = url }
        else { images.append(url); astro_add_image(s, url.path) }
    }
    func clearImages() { images = []; astro_clear_images(s); result = nil; known = []; identified = [:] }


    // ---------------------------------------------------------------- catalog
    func downloadCatalog() async {
        // center: the typed one; otherwise the plate center of the last reduction
        var r = ra.trimmingCharacters(in: .whitespaces), d = dec.trimmingCharacters(in: .whitespaces)
        if r.isEmpty || d.isEmpty {
            var cra = 0.0, cdec = 0.0
            if astro_frame_center(s, 0, &cra, &cdec) != 0 { r = Self.fmt(cra / 15, sign: false); d = Self.fmt(cdec, sign: true); ra = r; dec = d }
        }
        guard !r.isEmpty, !d.isEmpty else { status = L("Enter the field center (RA and Dec) to download Gaia DR3."); return }
        let base = String(cString: astro_catalog_url(s, r, d))
        guard !base.isEmpty else { status = L("Invalid center."); return }
        busy = true; defer { busy = false }
        for host in Online.vizierHosts {
            if Task.isCancelled { return }
            status = L("Downloading Gaia DR3 from %@…", host)
            guard let url = URL(string: base.replacingOccurrences(of: "vizier.cds.unistra.fr", with: host)) else { continue }
            guard let data = try? await Online.get(url, minBytes: 200, tries: 2, timeout: 40) else { continue }
            let text = String(decoding: data, as: UTF8.self)
            let rows = text.split(separator: "\n").filter { $0.first.map { $0.isNumber || $0 == " " } ?? false }.count
            if rows < 5 { continue }                                   // empty/error response: try the next mirror
            let dst = FileManager.default.temporaryDirectory.appendingPathComponent("gaia.tsv")
            do { try data.write(to: dst) } catch { status = L("Could not save the catalog: %@", error.localizedDescription); return }
            catalog = dst
            status = L("Gaia DR3 catalog downloaded from %@ (%d kB).", host, data.count / 1024); return
        }
        status = L("No VizieR mirror answered (no internet?). Try again in a moment.")
    }

    // ---------------------------------------------------------------- long-running tasks (cancellable)
    func start(_ op: @escaping () async -> Void) {
        work?.cancel()
        work = Task { await op() }
    }
    func cancel() {
        work?.cancel(); work = nil; busy = false; status = L("Cancelled.")
    }


    static func fmt(_ v: Double, sign: Bool) -> String {
        // Round to tenths of a second before splitting (never "59 60.0"); RA (unsigned) modulo 24.
        guard v.isFinite else { return sign ? "+00 00 00.0" : "00 00 00.0" }
        var t = Int((abs(v) * 36000).rounded())
        var neg = v < 0
        if !sign { let day = 24 * 36000; t = neg ? (day - t % day) % day : t % day; neg = false }
        if t == 0 { neg = false }
        let h = t / 36000, m = (t / 600) % 60, ds = t % 600
        return (sign ? (neg ? "-" : "+") : "") + String(format: "%02d %02d %02d.%d", h, m, ds / 10, ds % 10)
    }

    // ---------------------------------------------------------------- automatic pipeline
    private var sessionID = 0
    /// Back to home: switch screens first, then (after the animation) clear the data — the outgoing screen is still drawn.
    func newSession() {
        work?.cancel(); work = nil; busy = false; needCenter = false
        inspectPoint = nil; inspection = nil
        sessionID += 1; let id = sessionID
        phase = .home
        DispatchQueue.main.asyncAfter(deadline: .now() + 0.6) { [weak self] in
            guard let self, self.sessionID == id, self.phase == .home else { return }
            self.clearImages(); self.catalog = nil; self.ra = ""; self.dec = ""; self.object = ""; if self.isExample { self.restoreDefaults() }; self.isExample = false; self.steps = Self.freshSteps()
        }
    }

    func startImport(_ urls: [URL]) {
        sessionID += 1
        clearImages(); catalog = nil; ra = ""; dec = ""; object = ""; if isExample { restoreDefaults() }; isExample = false
        for u in urls { add(u) }
        applyHeader(); steps = Self.freshSteps(); phase = .running
        start { await self.runPipeline() }
    }

    /// Built-in example (Example/ in the bundle): example.cfg provides settings, images, catalog and center; runs offline.
    func openExample() {
        guard let dir = Bundle.main.url(forResource: "Example", withExtension: nil) else { status = L("Example data not found in the app bundle."); return }
        sessionID += 1
        clearImages(); catalog = nil; object = ""
        let cfg = dir.appendingPathComponent("example.cfg")
        astro_load_cfg(s, cfg.path); cfgName = "example.cfg"
        for n in get("example_images").split(separator: ";") {
            let u = dir.appendingPathComponent(n.trimmingCharacters(in: .whitespaces)); images.append(u); astro_add_image(s, u.path)
        }
        catalog = dir.appendingPathComponent(get("example_catalog"))
        ra = get("example_ra"); dec = get("example_dec")
        isExample = true; profileName = L("Example (offline)")
        steps = Self.freshSteps(); phase = .running
        start { await self.runPipeline() }
    }

    /// Reads the center (OBJCTRA/OBJCTDEC, RA/DEC, CRVAL1/2) and equipment (FOCALLEN, XPIXSZ/YPIXSZ) from the FITS header.
    func applyHeader() {
        guard let u = images.first, let h = FITSHeader.read(u) else { return }
        if ra.isEmpty || dec.isEmpty {
            if let r = h["OBJCTRA"], let d = h["OBJCTDEC"] { ra = r; dec = d }
            else if let r = h["RA"].flatMap(Double.init), let d = h["DEC"].flatMap(Double.init) { ra = Self.fmt(r / 15, sign: false); dec = Self.fmt(d, sign: true) }
            else if let r = h["CRVAL1"].flatMap(Double.init), let d = h["CRVAL2"].flatMap(Double.init) { ra = Self.fmt(r / 15, sign: false); dec = Self.fmt(d, sign: true) }
        }
        if let f = h["FOCALLEN"].flatMap(Double.init), f > 0 { set("focal_length_mm", String(f)) }
        if let x = h["XPIXSZ"].flatMap(Double.init), x > 0 { set("pixel_width_mm", String(x / 1000)); set("pixel_height_mm", String((h["YPIXSZ"].flatMap(Double.init) ?? x) / 1000)) }
    }

    /// Resolves a name (M42, NGC 2264, Betelgeuse...) via CDS Sesame.
    func resolve(_ name: String) async -> Bool {
        let q = name.addingPercentEncoding(withAllowedCharacters: .alphanumerics) ?? name
        for base in ["https://vizier.cfa.harvard.edu/viz-bin/nph-sesame/-oI/A?", "https://cds.unistra.fr/cgi-bin/nph-sesame/-oI/A?"] {
            guard let url = URL(string: base + q), let d = try? await Online.get(url, tries: 1, timeout: 20) else { continue }
            for l in String(decoding: d, as: UTF8.self).split(separator: "\n") where l.hasPrefix("%J ") {
                let v = l.dropFirst(3).split(separator: " ").compactMap { Double($0) }
                if v.count >= 2 { ra = Self.fmt(v[0] / 15, sign: false); dec = Self.fmt(v[1], sign: true); return true }
            }
        }
        return false
    }

    func runPipeline() async {
        phase = .running
        guard !images.isEmpty else { step(0, .failed(coreError(code: "no_images", message: nil, detail: nil))); return }
        step(0, .done(images.count >= 3 ? L("%d images", images.count) : L("%d image(s) — motion needs 3+", images.count)))
        if ra.isEmpty || dec.isEmpty { applyHeader() }
        if ra.isEmpty || dec.isEmpty { step(1, .failed(L("Enter the field center"))); needCenter = true; return }
        if catalog == nil {
            step(1, .running); await downloadCatalog()
            guard catalog != nil else { step(1, .failed(status)); return }
        }
        step(1, .done(catalog?.lastPathComponent ?? "")); if Task.isCancelled { return }
        if Task.isCancelled || phase != .running { return }
        step(2, .running); try? await Task.sleep(nanoseconds: 150_000_000)
        if Task.isCancelled || phase != .running { return }
        reduce()
        guard let fs = result?.frames, !fs.isEmpty, result?.error == nil else { step(2, .failed(result?.errorText ?? status)); return }
        let rms = fs.map { (($0.dra ?? 0) + ($0.dde ?? 0)) / 2 }.reduce(0, +) / Double(fs.count)
        step(2, .done(L("%d refs per image · %.2f″", fs.map(\.nref).reduce(0, +) / fs.count, rms)))
        step(3, .done(movers.isEmpty ? L("none") : L("%d object(s)", movers.count)))
        if Task.isCancelled || phase != .running { return }
        step(4, .running); await identify()
        step(4, .done(movers.isEmpty ? L("%d known asteroids in the field", known.count) : (movers.compactMap { identified[$0.desig] }.first.map { String($0.split(separator: "·").first ?? "") } ?? L("not identified"))))
        if !Task.isCancelled && phase == .running { phase = .results }
    }

    // ---------------------------------------------------------------- reduction and identification
    func reduce() {
        guard let cat = catalog else { status = L("Catalog missing (download Gaia DR3 or import one)."); return }
        busy = true
        let t0 = Date()
        let n = astro_reduce(s, cat.path, ra, dec, object)
        let json = String(cString: astro_result_json(s))
        result = try? JSONDecoder().decode(ReduceResult.self, from: Data(json.utf8))
        busy = false; known = []; identified = [:]
        if n < 0 { status = result?.errorText ?? L("Reduction error."); return }
        status = L("Reduction finished in %.2f s: %d images, %d moving object(s).", Date().timeIntervalSince(t0), frames.count, Int(n))
        if phase != .running { Task { await identify() } }
    }

    /// Identifies each moving object and lists known asteroids in the field (IMCCE SkyBoT, topocentric via the MPC code).
    func identify() async {
        guard let f0 = frames.first else { return }
        let code = get("mpc_code")
        var cra = 0.0, cdec = 0.0
        if astro_frame_center(s, 0, &cra, &cdec) != 0 {
            let radius = Double(max(f0.w, f0.h)) * 1.1 / 3600 * 0.75
            if let list = try? await Online.skybot(jd: f0.jd, ra: cra, dec: cdec, radiusDeg: radius, observer: code) { known = list }
        }
        guard var ms = result?.movers else { status += L(" Field: %d known asteroid(s).", known.count); return }
        for i in ms.indices {
            let q = ms[i].meas[0]
            guard let list = try? await Online.skybot(jd: q.jd, ra: q.rad, dec: q.decd, radiusDeg: 2.0 / 60, observer: code), let a = list.first else {
                identified[ms[i].desig] = L("not identified (new object or SkyBoT unavailable)"); continue
            }
            let desig = object.isEmpty ? (a.number ?? a.name) : object
            let label = (a.number.map { "(\($0)) " } ?? "") + a.name
            astro_set_mover_desig(s, Int32(i), desig)
            identified[desig] = String(format: "%@ · %@ · V %.1f · %.0f″", label, a.cls, a.vmag, a.sep)
            ms[i].desig = desig
        }
        result?.movers = ms
        status = L("Identification finished. %d known asteroid(s) in the field.", known.count)
    }

    @Published var inspection: Online.Inspection?
    @Published var inspectPoint: (x: Double, y: Double, ra: Double, dec: Double, frame: Int)?
    @Published var inspecting = false
    func inspect(frame k: Int, x: Double, y: Double) {
        var r = 0.0, d = 0.0
        guard astro_pix_to_sky(s, Int32(k), x, y, &r, &d) != 0, k < frames.count else { return }
        inspectPoint = (x, y, r, d, k); inspection = nil; inspecting = true
        UISelectionFeedbackGenerator().selectionChanged()
        let jd = frames[k].jd, code = get("mpc_code")
        Task { let res = await Online.whatIsHere(ra: r, dec: d, jd: jd, observer: code); if self.inspectPoint?.x == x { self.inspection = res; self.inspecting = false } }
    }

    func pix(frame k: Int, ra: Double, dec: Double) -> CGPoint? {
        var x = 0.0, y = 0.0
        return astro_sky_to_pix(s, Int32(k), ra, dec, &x, &y) != 0 ? CGPoint(x: x, y: y) : nil
    }

    func image(frame k: Int) -> CGImage? {
        var w: Int32 = 0, h: Int32 = 0
        guard astro_frame_size(s, Int32(k), &w, &h) != 0 else { return nil }
        var buf = [UInt8](repeating: 0, count: Int(w) * Int(h) * 4)
        astro_render(s, Int32(k), &buf)
        let prov = CGDataProvider(data: Data(buf) as CFData)!
        return CGImage(width: Int(w), height: Int(h), bitsPerComponent: 8, bitsPerPixel: 32, bytesPerRow: Int(w) * 4,
                       space: CGColorSpaceCreateDeviceRGB(), bitmapInfo: CGBitmapInfo(rawValue: CGImageAlphaInfo.noneSkipLast.rawValue),
                       provider: prov, decode: nil, shouldInterpolate: false, intent: .defaultIntent)
    }

    func reportFiles() -> [URL] {
        let dir = FileManager.default.temporaryDirectory
        return [("ADESReport.psv", 0), ("PhotReport.txt", 1), ("MPCReport.txt", 2)].map { name, k in
            let u = dir.appendingPathComponent(name)
            try? String(cString: astro_report(s, Int32(k))).write(to: u, atomically: true, encoding: .utf8)
            return u
        }
    }
}

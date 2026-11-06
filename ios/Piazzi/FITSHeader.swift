import Foundation

/// Minimal FITS header reader (80-character cards up to END).
enum FITSHeader {
    static func read(_ url: URL) -> [String: String]? {
        guard let fh = try? FileHandle(forReadingFrom: url) else { return nil }
        defer { try? fh.close() }
        guard let data = try? fh.read(upToCount: 2880 * 12), data.count >= 2880,
              String(decoding: data.prefix(6), as: UTF8.self) == "SIMPLE" else { return nil }
        let text = String(decoding: data, as: UTF8.self); var out: [String: String] = [:]
        var i = text.startIndex
        while i < text.endIndex, let j = text.index(i, offsetBy: 80, limitedBy: text.endIndex) {
            let card = String(text[i..<j]); i = j
            let key = card.prefix(8).trimmingCharacters(in: .whitespaces)
            if key == "END" { break }
            guard card.count > 10, card.dropFirst(8).hasPrefix("= ") else { continue }
            var v = String(card.dropFirst(10))
            if let q = v.firstIndex(of: "'") {
                let rest = v[v.index(after: q)...]; v = String(rest.prefix { $0 != "'" })
            } else if let sl = v.firstIndex(of: "/") { v = String(v[..<sl]) }
            out[key] = v.trimmingCharacters(in: .whitespaces)
        }
        return out
    }
}

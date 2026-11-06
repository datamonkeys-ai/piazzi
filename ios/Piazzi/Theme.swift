import SwiftUI

/// "Pocket observatory" identity: black ink, instrument amber (preserves night vision),
/// serif titles, monospaced data, star chart background.
enum Theme {
    static let ink = Color(red: 0.027, green: 0.031, blue: 0.039)
    static let panel = Color(red: 0.055, green: 0.063, blue: 0.075)
    static let line = Color.white.opacity(0.09)
    static let amber = Color(red: 1.0, green: 0.71, blue: 0.28)
    static let amberDim = Color(red: 1.0, green: 0.71, blue: 0.28).opacity(0.16)
    static let red = Color(red: 1.0, green: 0.36, blue: 0.31)
    static let cyan = Color(red: 0.50, green: 0.83, blue: 1.0)
    static let green = Color(red: 0.61, green: 0.90, blue: 0.39)
    static let text = Color(red: 0.93, green: 0.91, blue: 0.86)
    static let dim = Color(red: 0.93, green: 0.91, blue: 0.86).opacity(0.55)
    static let faint = Color(red: 0.93, green: 0.91, blue: 0.86).opacity(0.32)

    // aliases used by other screens
    static var gold: Color { amber }
    static var good: Color { green }
    static var warn: Color { amber }
    static var bad: Color { red }
    static var sky: Color { cyan }
    static var card: Color { panel }

    static func title(_ size: CGFloat) -> Font { .system(size: size, weight: .semibold, design: .serif) }
    static func mono(_ size: CGFloat, _ w: Font.Weight = .regular) -> Font { .system(size: size, weight: w, design: .monospaced) }
    static var backdrop: some View { ZStack { ink; StarChart() }.ignoresSafeArea() }
}

/// Instrument label: uppercase, monospaced, letter-spaced.
struct Tag: View {
    let text: String; var color: Color = Theme.dim
    init(_ t: String, color: Color = Theme.dim) { text = t; self.color = color }
    var body: some View { Text(text.uppercased()).font(Theme.mono(10.5, .medium)).tracking(1.6).foregroundStyle(color) }
}

/// Background: star chart (curved RA/Dec grid + stars + vignette).
struct StarChart: View {
    var body: some View {
        Canvas { ctx, size in
            let w = size.width, h = size.height
            var grid = Path()
            for i in 0...8 { let x = w * CGFloat(i) / 8; grid.move(to: CGPoint(x: x, y: 0)); grid.addQuadCurve(to: CGPoint(x: x + (x - w / 2) * 0.12, y: h), control: CGPoint(x: x + (x - w / 2) * 0.03, y: h / 2)) }
            for j in 0...12 { let y = h * CGFloat(j) / 12; grid.move(to: CGPoint(x: 0, y: y)); grid.addQuadCurve(to: CGPoint(x: w, y: y), control: CGPoint(x: w / 2, y: y - 18)) }
            ctx.stroke(grid, with: .color(.white.opacity(0.035)), lineWidth: 0.6)
            var rng = Lcg(seed: 7)
            for _ in 0..<220 {
                let x = CGFloat(rng.next() % 10_000) / 10_000 * w, y = CGFloat(rng.next() % 10_000) / 10_000 * h
                let m = CGFloat(rng.next() % 1000) / 1000, r = 0.35 + pow(m, 6) * 2.4
                ctx.fill(Path(ellipseIn: CGRect(x: x, y: y, width: r, height: r)), with: .color(.white.opacity(0.15 + Double(m) * 0.5)))
            }
            ctx.fill(Path(CGRect(origin: .zero, size: size)), with: .radialGradient(Gradient(colors: [.clear, .black.opacity(0.65)]),
                center: CGPoint(x: w / 2, y: h * 0.35), startRadius: min(w, h) * 0.25, endRadius: max(w, h) * 0.8))
        }.allowsHitTesting(false)
    }
}
struct Lcg { var s: UInt64; init(seed: UInt64) { s = seed }; mutating func next() -> UInt64 { s = s &* 6364136223846793005 &+ 1442695040888963407; return s >> 33 } }

/// Reticle corner marks around a point.
func reticle(_ c: CGPoint, _ r: CGFloat, gap: CGFloat = 0.45) -> Path {
    var p = Path(); let l = r * gap
    for (sx, sy) in [(-1.0, -1.0), (1, -1), (-1, 1), (1, 1)] {
        let cx = c.x + r * sx, cy = c.y + r * sy
        p.move(to: CGPoint(x: cx, y: cy - l * sy)); p.addLine(to: CGPoint(x: cx, y: cy)); p.addLine(to: CGPoint(x: cx - l * sx, y: cy))
    }
    return p
}

/// Primary amber button.
struct PrimaryButtonStyle: ButtonStyle {
    func makeBody(configuration: Configuration) -> some View {
        configuration.label.font(Theme.mono(14, .bold)).tracking(1.4).textCase(.uppercase).foregroundStyle(Theme.ink)
            .frame(maxWidth: .infinity).padding(.vertical, 15)
            .background(Theme.amber.opacity(configuration.isPressed ? 0.75 : 1), in: RoundedRectangle(cornerRadius: 6))
    }
}
struct GhostButtonStyle: ButtonStyle {
    func makeBody(configuration: Configuration) -> some View {
        configuration.label.font(Theme.mono(12.5, .semibold)).tracking(1.2).textCase(.uppercase).foregroundStyle(Theme.text)
            .frame(maxWidth: .infinity).padding(.vertical, 13)
            .overlay(RoundedRectangle(cornerRadius: 6).stroke(Theme.line, lineWidth: 1))
            .opacity(configuration.isPressed ? 0.6 : 1)
    }
}

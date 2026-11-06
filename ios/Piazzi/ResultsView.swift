import Combine
import SwiftUI

struct ResultsView: View {
    @EnvironmentObject var e: Engine
    @Environment(\.horizontalSizeClass) private var hsize
    @State private var k = 0
    @State private var blink = true
    @State private var showSheet = true
    @State private var detent: PresentationDetent = .height(230)
    @State private var layers = Layers()
    @State private var zoom: CGFloat = 1
    @GestureState private var pinch: CGFloat = 1
    let timer = Timer.publish(every: 0.6, on: .main, in: .common).autoconnect()
    struct Layers { var refs = false; var dets = false; var known = true; var track = true }

    var body: some View {
        if hsize == .regular { wide } else { compact }
    }

    /// Mac / iPad: image and card side by side.
    var wide: some View {
        HStack(spacing: 0) {
            ZStack(alignment: .top) {
                Theme.ink.ignoresSafeArea()
                if !e.frames.isEmpty {
                    SkyImage(k: k, layers: layers, scale: zoom * pinch) { x, y in e.inspect(frame: k, x: x, y: y); blink = false }
                        .gesture(MagnificationGesture().updating($pinch) { v, st, _ in st = v }.onEnded { zoom = max(1, min(6, zoom * $0)) })
                }
                topBar
            }
            Rectangle().fill(Theme.line).frame(width: 1)
            ResultsSheet(k: $k).frame(width: 420).background(Theme.ink)
        }
        .onReceive(timer) { _ in if blink, !e.frames.isEmpty, e.phase == .results { k = (k + 1) % e.frames.count } }
        .onAppear { showSheet = false }
    }

    var compact: some View {
        ZStack(alignment: .top) {
            Theme.ink.ignoresSafeArea()
            if !e.frames.isEmpty {
                SkyImage(k: k, layers: layers, scale: zoom * pinch) { x, y in
                    e.inspect(frame: k, x: x, y: y); blink = false
                    if detent == .height(230) { detent = .medium }
                }
                .gesture(MagnificationGesture().updating($pinch) { v, st, _ in st = v }.onEnded { zoom = max(1, min(6, zoom * $0)) })
            }
            topBar
        }
        .onReceive(timer) { _ in if blink, !e.frames.isEmpty, e.phase == .results { k = (k + 1) % e.frames.count } }
        .sheet(isPresented: $showSheet) {
            ResultsSheet(k: $k)
                .presentationDetents([.height(230), .medium, .large], selection: $detent)
                .presentationBackgroundInteraction(.enabled(upThrough: .medium))
                .presentationBackground(Theme.ink.opacity(0.97))
                .presentationCornerRadius(14)
                .presentationDragIndicator(.visible)
                .interactiveDismissDisabled()
                .preferredColorScheme(.dark)
        }
        .onAppear { UINotificationFeedbackGenerator().notificationOccurred(e.movers.isEmpty ? .warning : .success) }
    }

    func chip<L: View>(@ViewBuilder _ l: () -> L) -> some View {
        l().frame(width: 38, height: 38).background(Theme.panel.opacity(0.85), in: RoundedRectangle(cornerRadius: 6)).overlay(RoundedRectangle(cornerRadius: 6).stroke(Theme.line))
    }

    var topBar: some View {
        VStack(spacing: 8) {
            HStack(spacing: 10) {
                Button { showSheet = false; blink = false; e.newSession() } label: { chip { Image(systemName: "chevron.left").font(.system(size: 15, weight: .semibold)) } }
                Spacer()
                HStack(spacing: 4) {
                    ForEach(0..<e.frames.count, id: \.self) { i in
                        Button { blink = false; k = i } label: {
                            Text("\(i + 1)").font(Theme.mono(12, .bold)).frame(width: 26, height: 26)
                                .foregroundStyle(i == k ? Theme.ink : Theme.dim).background(i == k ? Theme.amber : .clear, in: RoundedRectangle(cornerRadius: 4))
                        }
                    }
                }.padding(4).background(Theme.panel.opacity(0.85), in: RoundedRectangle(cornerRadius: 6)).overlay(RoundedRectangle(cornerRadius: 6).stroke(Theme.line))
                Spacer()
                Button { blink.toggle() } label: { chip { Image(systemName: blink ? "pause.fill" : "play.fill").font(.system(size: 13)) } }
                Menu {
                    Toggle(L("Object trail"), isOn: $layers.track)
                    Toggle(L("Known asteroids"), isOn: $layers.known)
                    Toggle(L("Reference stars"), isOn: $layers.refs)
                    Toggle(L("All detections"), isOn: $layers.dets)
                } label: { chip { Image(systemName: "square.3.layers.3d").font(.system(size: 14)) } }
            }
            if k < e.frames.count {
                HStack {
                    Tag(L("Plate %d · %@", k + 1, String(e.frames[k].time.suffix(8))), color: Theme.amber).lineLimit(1)
                    Spacer()
                    if let p = e.inspectPoint { Tag("\(Engine.fmt(p.ra / 15, sign: false))  \(Engine.fmt(p.dec, sign: true))", color: Theme.dim) }
                    else { Tag(L("tap = identify"), color: Theme.faint).lineLimit(1) }
                }.padding(.horizontal, 4)
            }
        }
        .foregroundStyle(Theme.text).padding(.horizontal, 14).padding(.top, 6)
    }
}

/// Image with overlays (reticles) drawn on top.
struct SkyImage: View {
    @EnvironmentObject var e: Engine
    let k: Int; let layers: ResultsView.Layers; let scale: CGFloat; let onTap: (Double, Double) -> Void
    var body: some View {
        GeometryReader { geo in
          if e.frames.isEmpty { Color.clear } else {
            let f = e.frames[max(0, min(k, e.frames.count - 1))]
            let k = max(0, min(k, e.frames.count - 1))
            let base = min(geo.size.width / CGFloat(f.w), (geo.size.height * 0.62) / CGFloat(f.h))
            let s = base * scale
            ScrollView([.horizontal, .vertical], showsIndicators: false) {
                ZStack(alignment: .topLeading) {
                    if let cg = e.image(frame: k) {
                        Image(decorative: cg, scale: 1).interpolation(.none).resizable().frame(width: CGFloat(f.w) * s, height: CGFloat(f.h) * s)
                    }
                    Canvas { ctx, _ in draw(ctx, f, s, k) }.frame(width: CGFloat(f.w) * s, height: CGFloat(f.h) * s)
                }
                .contentShape(Rectangle())
                .onTapGesture { loc in
                    var x = Double(loc.x / s) - 0.5, y = Double(loc.y / s) - 0.5
                    if let d = f.dets.min(by: { hypot($0[0] - x, $0[1] - y) < hypot($1[0] - x, $1[1] - y) }), hypot(d[0] - x, d[1] - y) < 8 { x = d[0]; y = d[1] }
                    onTap(x, y)
                }
                .frame(minWidth: geo.size.width, minHeight: geo.size.height * 0.62, alignment: .top)
            }
            .defaultScrollAnchor(.top)
            .padding(.top, 84)
          }
        }
    }

    func draw(_ ctx: GraphicsContext, _ f: FrameResult, _ s: CGFloat, _ k: Int) {
        func p(_ x: Double, _ y: Double) -> CGPoint { CGPoint(x: (x + 0.5) * s, y: (y + 0.5) * s) }
        if layers.dets { for d in f.dets { let c = p(d[0], d[1]); ctx.stroke(Path(ellipseIn: CGRect(x: c.x - 3 * s, y: c.y - 3 * s, width: 6 * s, height: 6 * s)), with: .color(Theme.cyan.opacity(0.35)), lineWidth: 0.8) } }
        if layers.refs { for r in f.refs { ctx.stroke(reticle(p(r[0], r[1]), 5 * s, gap: 0.5), with: .color(r[2] > 0 ? Theme.green.opacity(0.7) : Theme.red.opacity(0.7)), lineWidth: 1) } }
        if layers.known {
            for a in e.known where !e.movers.contains(where: { e.identified[$0.desig]?.contains(a.name) == true }) {
                guard let q = e.pix(frame: k, ra: a.ra, dec: a.dec) else { continue }
                let c = p(q.x, q.y)
                ctx.stroke(reticle(c, 7 * s, gap: 0.35), with: .color(Theme.cyan.opacity(0.7)), lineWidth: 1)
                ctx.draw(Text(a.name).font(Theme.mono(9, .medium)).foregroundColor(Theme.cyan.opacity(0.85)), at: CGPoint(x: c.x, y: c.y - 11 * s))
            }
        }
        for m in e.movers {
            let pts = m.meas.map { p($0.x, $0.y) }
            if layers.track, pts.count > 1 {
                var path = Path(); path.move(to: pts[0]); for q in pts.dropFirst() { path.addLine(to: q) }
                ctx.stroke(path, with: .color(Theme.amber.opacity(0.5)), style: StrokeStyle(lineWidth: 1.2, dash: [3, 3]))
                for q in pts { ctx.fill(Path(ellipseIn: CGRect(x: q.x - 1.5, y: q.y - 1.5, width: 3, height: 3)), with: .color(Theme.amber.opacity(0.6))) }
            }
            if k < pts.count {
                ctx.stroke(reticle(pts[k], 10 * s), with: .color(Theme.amber), lineWidth: 1.8)
                let label = e.identified[m.desig].map { String($0.split(separator: "·").first ?? "") } ?? m.label
                ctx.draw(Text(label.trimmingCharacters(in: .whitespaces)).font(Theme.mono(11, .bold)).foregroundColor(Theme.amber), at: CGPoint(x: pts[k].x, y: pts[k].y + 18 * s))
            }
        }
        if let ip = e.inspectPoint { ctx.stroke(reticle(p(ip.x, ip.y), 8 * s, gap: 0.6), with: .color(.white), lineWidth: 1.4) }
    }
}

// MARK: - panel
struct ResultsSheet: View {
    @EnvironmentObject var e: Engine
    @Binding var k: Int
    var body: some View {
        ScrollView(showsIndicators: false) {
            VStack(alignment: .leading, spacing: 14) {
                if e.inspectPoint != nil { InspectCard() }
                if e.movers.isEmpty {
                    VStack(alignment: .leading, spacing: 6) {
                        Tag(L("Result"), color: Theme.amber)
                        Text(L("Nothing moved.")).font(Theme.title(26))
                        Text(L("The plates were solved, but no object moved between the images above the threshold. %d known asteroid(s) cross the field, probably too faint.", e.known.count)).font(.system(size: 14)).foregroundStyle(Theme.dim)
                    }
                }
                ForEach(e.movers) { m in ObjectCard(m: m) }
                QualityCard()
                if !e.known.isEmpty { KnownCard() }
                Button(L("New session")) { e.newSession() }.buttonStyle(GhostButtonStyle())
            }
            .padding(.horizontal, 20).padding(.top, 22).padding(.bottom, 30)
        }
    }
}

struct Row: View {
    let k: String, v: String; var color: Color = Theme.text
    var body: some View {
        HStack(alignment: .firstTextBaseline) { Tag(k); Spacer(); Text(v).font(Theme.mono(13.5, .medium)).foregroundStyle(color).multilineTextAlignment(.trailing) }
            .padding(.vertical, 7).overlay(alignment: .bottom) { Rectangle().fill(Theme.line).frame(height: 1) }
    }
}

struct Panel<C: View>: View {
    let content: C; var accent: Color? = nil
    init(accent: Color? = nil, @ViewBuilder _ c: () -> C) { content = c(); self.accent = accent }
    var body: some View {
        VStack(alignment: .leading, spacing: 10) { content }.padding(16)
            .background(Theme.panel, in: RoundedRectangle(cornerRadius: 8))
            .overlay(RoundedRectangle(cornerRadius: 8).stroke(accent?.opacity(0.45) ?? Theme.line, lineWidth: 1))
    }
}

struct ObjectCard: View {
    @EnvironmentObject var e: Engine
    let m: Mover
    var id: String? { e.identified[m.desig] }
    var known: Bool { id.map { !$0.hasPrefix(L("not identified")) } ?? false }
    var parts: [String] { id?.split(separator: "·").map { $0.trimmingCharacters(in: .whitespaces) } ?? [] }
    var body: some View {
        Panel(accent: known ? Theme.amber : Theme.red) {
            HStack { Tag(known ? L("Object card · identified") : L("Object card · not identified"), color: known ? Theme.amber : Theme.red); Spacer()
                Image(systemName: known ? "checkmark.seal.fill" : "questionmark.diamond.fill").foregroundStyle(known ? Theme.amber : Theme.red) }
            Text(known ? (parts.first ?? m.desig) : L("Possible new object")).font(Theme.title(30)).lineLimit(2).minimumScaleFactor(0.7)
            VStack(spacing: 0) {
                if known, parts.count > 1 { Row(k: L("Class"), v: parts[1]) }
                if known, parts.count > 2 { Row(k: L("Predicted brightness"), v: parts[2]) }
                if known, parts.count > 3 { Row(k: L("Offset from ephemeris"), v: parts[3]) }
                if m.desig.isEmpty { Row(k: L("Temporary designation"), v: m.label) }
                Row(k: L("Rate"), v: String(format: "%.2f″/min", m.rate ?? 0))
                Row(k: L("Position angle"), v: String(format: "%.1f°", m.pa ?? 0))
                Row(k: L("Measured magnitude"), v: String(format: "%.1f", m.meas.map(\.mag).reduce(0, +) / Double(max(m.meas.count, 1))))
            }
            Tag(L("Measurements")).padding(.top, 6)
            VStack(spacing: 4) {
                ForEach(Array(m.meas.enumerated()), id: \.offset) { i, q in
                    HStack(spacing: 8) {
                        Text("\(i + 1)").font(Theme.mono(11, .bold)).foregroundStyle(Theme.amber).frame(width: 14)
                        Text(String(q.time.suffix(8))).font(Theme.mono(11.5)).foregroundStyle(Theme.dim)
                        Spacer()
                        Text("\(q.ra)  \(q.dec)").font(Theme.mono(11.5)).lineLimit(1).minimumScaleFactor(0.6)
                        Text(String(format: "%.1f", q.mag)).font(Theme.mono(11.5)).foregroundStyle(Theme.dim)
                    }
                }
            }
            if e.isExample {
                // built-in example: the data is not the user's — never submit to the MPC
                Button(L("Send to MPC")) {}.buttonStyle(PrimaryButtonStyle()).disabled(true).opacity(0.4).padding(.top, 8)
                Tag(L("Example data — not your observation"), color: Theme.red)
            } else {
                ShareLink(items: [e.reportFiles()[0]]) { Text(L("Send to MPC")) }.buttonStyle(PrimaryButtonStyle()).padding(.top, 8)
            }
            ShareLink(items: e.reportFiles()) { Text(L("ADES · photometric · MPC 80")) }.buttonStyle(GhostButtonStyle())
        }
    }
}

struct InspectCard: View {
    @EnvironmentObject var e: Engine
    var body: some View {
        Panel(accent: .white) {
            HStack { Tag(L("What is this?"), color: Theme.text); Spacer()
                Button { e.inspectPoint = nil; e.inspection = nil } label: { Image(systemName: "xmark").font(.system(size: 12, weight: .bold)) }.foregroundStyle(Theme.dim) }
            if let p = e.inspectPoint {
                Text("\(Engine.fmt(p.ra / 15, sign: false))   \(Engine.fmt(p.dec, sign: true))").font(Theme.mono(13)).foregroundStyle(Theme.dim)
            }
            if e.inspecting {
                HStack(spacing: 10) { ProgressView().tint(Theme.amber); Text(L("Querying SIMBAD, VSX, Gaia and SkyBoT…")).font(Theme.mono(12)).foregroundStyle(Theme.dim) }.padding(.vertical, 6)
            } else if let r = e.inspection {
                let v = r.verdict
                Text(v.title).font(Theme.title(24))
                Text(v.detail).font(.system(size: 14)).foregroundStyle(color(v.kind))
                VStack(spacing: 0) {
                    if let a = r.asteroid { Row(k: "SkyBoT", v: String(format: "%@ · %.0f″", a.name, a.sep)) }
                    ForEach(Array(r.vsx.prefix(2).enumerated()), id: \.offset) { _, x in Row(k: "VSX", v: "\(x.name) · \(x.type)") }
                    ForEach(Array(r.simbad.prefix(2).enumerated()), id: \.offset) { _, x in Row(k: "SIMBAD", v: String(format: "%@ · %@ · %.1f″", x.id, x.type, x.sep)) }
                    if let g = r.gaia { Row(k: "Gaia DR3", v: String(format: "G %.2f%@ · %.1f″", g.g, g.bprp.map { String(format: " · BP−RP %.2f", $0) } ?? "", g.sep)) }
                }
                if !r.failed.isEmpty { Tag(L("no answer: ") + r.failed.joined(separator: ", "), color: Theme.faint) }
            }
        }
    }
    func color(_ k: Online.Inspection.Kind) -> Color { switch k { case .asteroid: return Theme.amber; case .variable: return Theme.red; case .star: return Theme.green; case .other: return Theme.cyan; case .unknown: return Theme.dim } }
}

struct QualityCard: View {
    @EnvironmentObject var e: Engine
    var body: some View {
        Panel {
            Tag(L("Plate quality"))
            ForEach(Array(e.frames.enumerated()), id: \.offset) { i, f in
                let rms = ((f.dra ?? 0) + (f.dde ?? 0)) / 2
                let (txt, col): (String, Color) = f.nref < 6 ? (L("few stars"), Theme.red) : rms < 0.25 ? (L("excellent"), Theme.green) : rms < 0.6 ? (L("good"), Theme.amber) : (L("poor"), Theme.red)
                VStack(alignment: .leading, spacing: 5) {
                    HStack {
                        Text(L("Plate %d", i + 1)).font(Theme.mono(12.5, .semibold))
                        Text(txt).font(Theme.mono(12.5)).foregroundStyle(col)
                        Spacer()
                        Text(String(format: "%d refs · %.2f″ · %.1f″", f.nref, rms, f.seeing ?? 0)).font(Theme.mono(11)).foregroundStyle(Theme.dim)
                    }
                    GeometryReader { g in
                        ZStack(alignment: .leading) {
                            Rectangle().fill(Theme.line)
                            Rectangle().fill(col).frame(width: g.size.width * CGFloat(max(0.05, min(1, 1 - rms))))
                        }
                    }.frame(height: 3)
                }
            }
        }
    }
}

struct KnownCard: View {
    @EnvironmentObject var e: Engine
    var body: some View {
        Panel {
            Tag(L("Known asteroids in the field · SkyBoT"))
            VStack(spacing: 0) {
                ForEach(e.known.prefix(8)) { a in Row(k: (a.number.map { "(\($0)) " } ?? "") + a.name, v: String(format: "V %.1f", a.vmag), color: a.vmag < 18.5 ? Theme.cyan : Theme.dim) }
            }
            if e.known.count > 8 { Tag(L("+ %d fainter", e.known.count - 8), color: Theme.faint) }
        }
    }
}

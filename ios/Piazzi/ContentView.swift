import SwiftUI
import UniformTypeIdentifiers

// MARK: - root
struct ContentView: View {
    @EnvironmentObject var e: Engine
    @AppStorage("lang") private var lang = "system"
    var body: some View {
        ZStack {
            Theme.backdrop
            switch e.phase {
            case .home: HomeView().transition(.opacity)
            case .running: SessionView().transition(.opacity.combined(with: .offset(y: 12)))
            case .results: ResultsView().transition(.opacity)
            }
        }
        .animation(.easeInOut(duration: 0.4), value: e.phase)
        .preferredColorScheme(.dark)
        .tint(Theme.amber)
        .foregroundStyle(Theme.text)
        .sheet(isPresented: $e.needCenter) { CenterSheet() }
        .id(lang) // language switch: redraw everything with the new strings
    }
}

// MARK: - home
struct HomeView: View {
    @EnvironmentObject var e: Engine
    @State private var importing = false
    @State private var settings = false
    @State private var blink = false
    var body: some View {
        ScrollView(showsIndicators: false) {
            VStack(alignment: .leading, spacing: 0) {
                HStack {
                    Tag(L("Pocket observatory"))
                    Spacer()
                    Button { settings = true } label: {
                        HStack(spacing: 6) { Tag(e.profileName == "Default" ? L("Default") : e.profileName, color: Theme.amber); Image(systemName: "slider.horizontal.3").font(.system(size: 14)) }
                    }.foregroundStyle(Theme.amber)
                }
                .padding(.top, 6)

                Text("Piazzi").font(Theme.title(46)).padding(.top, 26)
                Text(L("Measure the sky. Find what moves.")).font(.system(size: 19, design: .serif)).italic().foregroundStyle(Theme.dim).padding(.top, 2)

                // sample plate with blinking reticle
                ZStack {
                    StarChart().clipShape(RoundedRectangle(cornerRadius: 4))
                    Canvas { ctx, size in
                        let c = CGPoint(x: size.width * (blink ? 0.58 : 0.55), y: size.height * (blink ? 0.43 : 0.47))
                        ctx.stroke(reticle(c, 15), with: .color(Theme.amber), lineWidth: 1.5)
                        ctx.fill(Path(ellipseIn: CGRect(x: c.x - 1.6, y: c.y - 1.6, width: 3.2, height: 3.2)), with: .color(.white))
                        var t = Path(); t.move(to: CGPoint(x: size.width * 0.49, y: size.height * 0.55)); t.addLine(to: CGPoint(x: size.width * 0.64, y: size.height * 0.35))
                        ctx.stroke(t, with: .color(Theme.amber.opacity(0.4)), style: StrokeStyle(lineWidth: 1, dash: [3, 4]))
                    }
                    VStack { Spacer(); HStack { Tag("07 11 44.28 +22 24 16", color: Theme.amber.opacity(0.8)).lineLimit(1); Spacer(); Tag("0.35″/min", color: Theme.amber.opacity(0.8)).lineLimit(1) } }.padding(10)
                }
                .frame(height: 150).overlay(RoundedRectangle(cornerRadius: 4).stroke(Theme.line))
                .padding(.top, 24)
                .onAppear { withAnimation(.easeInOut(duration: 0.7).repeatForever()) { blink.toggle() } }

                VStack(spacing: 12) {
                    PlateCard(n: "01", tag: L("Your plates"), title: L("Open my images"), detail: L("FITS or SBIG from your telescope. Catalog, plate, motion and identification — automatic."), primary: true) { importing = true }
                    PlateCard(n: "02", tag: L("Try it"), title: L("Open example"), detail: L("Three ZTF images of asteroid (5157) Hindemith with a Gaia DR3 catalog — runs fully offline.")) { e.openExample() }
                }.padding(.top, 26)

                HStack(spacing: 14) {
                    Tag("open source"); Tag("Gaia DR3"); Tag("SkyBoT")
                }.padding(.top, 26).foregroundStyle(Theme.faint)
                Text("Piazzi \(String(cString: astro_version())) · open source · MIT").font(Theme.mono(11)).foregroundStyle(Theme.faint).padding(.top, 10)
            }
            .padding(.horizontal, 22).padding(.bottom, 30)
        }
        .fileImporter(isPresented: $importing, allowedContentTypes: [.data, .item], allowsMultipleSelection: true) { r in
            guard case .success(let urls) = r else { return }
            var local: [URL] = []
            for u in urls {
                let ok = u.startAccessingSecurityScopedResource(); defer { if ok { u.stopAccessingSecurityScopedResource() } }
                let dst = FileManager.default.temporaryDirectory.appendingPathComponent(u.lastPathComponent)
                try? FileManager.default.removeItem(at: dst); try? FileManager.default.copyItem(at: u, to: dst); local.append(dst)
            }
            e.startImport(local)
        }
        .sheet(isPresented: $settings) { SettingsSheet() }
    }
}

struct PlateCard: View {
    let n: String, tag: String, title: String, detail: String; var primary = false; let action: () -> Void
    var body: some View {
        Button(action: { UIImpactFeedbackGenerator(style: .light).impactOccurred(); action() }) {
            HStack(alignment: .top, spacing: 16) {
                Text(n).font(Theme.mono(13, .bold)).foregroundStyle(primary ? Theme.ink : Theme.amber)
                    .frame(width: 34, height: 34).background(primary ? Theme.amber : Theme.amberDim, in: RoundedRectangle(cornerRadius: 4))
                VStack(alignment: .leading, spacing: 5) {
                    Tag(tag)
                    Text(title).font(Theme.title(21)).foregroundStyle(Theme.text)
                    Text(detail).font(.system(size: 14)).foregroundStyle(Theme.dim).multilineTextAlignment(.leading).fixedSize(horizontal: false, vertical: true)
                }
                Spacer(minLength: 0)
                Image(systemName: "arrow.up.right").font(.system(size: 14, weight: .semibold)).foregroundStyle(primary ? Theme.amber : Theme.faint).padding(.top, 2)
            }
            .padding(16)
            .background(Theme.panel.opacity(0.92), in: RoundedRectangle(cornerRadius: 8))
            .overlay(RoundedRectangle(cornerRadius: 8).stroke(primary ? Theme.amber.opacity(0.45) : Theme.line, lineWidth: 1))
        }.buttonStyle(.plain)
    }
}

// MARK: - session: observing log
struct SessionView: View {
    @EnvironmentObject var e: Engine
    @State private var t0 = Date()
    var failed: Bool { e.steps.contains { if case .failed = $0.state { return true }; return false } }
    var done: Int { e.steps.filter { if case .done = $0.state { return true }; return false }.count }
    var body: some View {
        VStack(alignment: .leading, spacing: 0) {
            HStack {
                Button { e.cancel(); e.newSession() } label: { Image(systemName: "xmark").font(.system(size: 15, weight: .semibold)) }.foregroundStyle(Theme.dim)
                Spacer()
                Tag(L("Step %d/%d", min(done + 1, e.steps.count), e.steps.count), color: Theme.amber)
            }.padding(.top, 8)
            // thin progress bar
            GeometryReader { g in
                ZStack(alignment: .leading) {
                    Rectangle().fill(Theme.line)
                    Rectangle().fill(failed ? Theme.red : Theme.amber).frame(width: g.size.width * CGFloat(done) / CGFloat(e.steps.count))
                }
            }.frame(height: 2).padding(.top, 14).animation(.easeInOut, value: done)

            Text(failed ? L("Needs attention.") : (done == e.steps.count ? L("Done.") : L("Reducing the plates…")))
                .font(Theme.title(32)).padding(.top, 28)
            Text(e.status).font(.system(size: 15)).foregroundStyle(Theme.dim).padding(.top, 6).lineLimit(3).contentTransition(.opacity)

            VStack(alignment: .leading, spacing: 0) {
                Tag(L("Observing log")).padding(.bottom, 12)
                ForEach(e.steps) { s in LogLine(step: s, t0: t0) }
            }
            .padding(16).background(Theme.panel.opacity(0.9), in: RoundedRectangle(cornerRadius: 8))
            .overlay(RoundedRectangle(cornerRadius: 8).stroke(Theme.line))
            .padding(.top, 26)

            Spacer()
            if failed {
                VStack(spacing: 10) {
                    Button(L("Try again")) { e.steps = Engine.freshSteps(); e.start { await e.runPipeline() } }.buttonStyle(PrimaryButtonStyle())
                    Button(L("Back to start")) { e.newSession() }.buttonStyle(GhostButtonStyle())
                }
            } else {
                Button(L("Cancel")) { e.cancel(); e.newSession() }.buttonStyle(GhostButtonStyle())
            }
        }
        .padding(.horizontal, 22).padding(.bottom, 16)
        .onAppear { t0 = Date() }
    }
}

struct LogLine: View {
    let step: Step; let t0: Date
    @State private var stamp = ""
    var body: some View {
        HStack(alignment: .firstTextBaseline, spacing: 10) {
            Text(stamp.isEmpty ? "--:--" : stamp).font(Theme.mono(11)).foregroundStyle(Theme.faint).frame(width: 40, alignment: .leading)
            glyph.frame(width: 16)
            VStack(alignment: .leading, spacing: 2) {
                Text(step.title).font(Theme.mono(13.5, .semibold)).foregroundStyle(step.state == .pending ? Theme.faint : Theme.text)
                if let d = detail { Text(d).font(Theme.mono(12)).foregroundStyle(isFail ? Theme.red : Theme.dim).fixedSize(horizontal: false, vertical: true) }
            }
            Spacer()
        }
        .padding(.vertical, 7)
        .onChange(of: step.state) { _, new in if new != .pending { let s = Int(Date().timeIntervalSince(t0)); stamp = String(format: "%02d:%02d", s / 60, s % 60) } }
        .animation(.easeOut(duration: 0.25), value: step.state)
    }
    var isFail: Bool { if case .failed = step.state { return true }; return false }
    @ViewBuilder var glyph: some View {
        switch step.state {
        case .pending: Text("○").font(Theme.mono(13)).foregroundStyle(Theme.faint)
        case .running: ProgressView().controlSize(.mini).tint(Theme.amber)
        case .done: Text("●").font(Theme.mono(13)).foregroundStyle(Theme.amber)
        case .failed: Text("✕").font(Theme.mono(13, .bold)).foregroundStyle(Theme.red)
        }
    }
    var detail: String? { switch step.state { case .done(let s), .failed(let s): return s.isEmpty ? nil : s; default: return nil } }
}

// MARK: - field center
struct CenterSheet: View {
    @EnvironmentObject var e: Engine
    @State private var name = ""
    @State private var searching = false
    @State private var msg = ""
    var body: some View {
        VStack(alignment: .leading, spacing: 18) {
            Tag(L("Field center"), color: Theme.amber)
            Text(L("Where was the telescope pointing?")).font(Theme.title(26))
            Text(L("I could not find the coordinates in the image header. Search by target name or type RA and Dec.")).font(.system(size: 14)).foregroundStyle(Theme.dim)
            HStack {
                TextField("M42, NGC 2264, Betelgeuse…", text: $name).font(Theme.mono(15)).textInputAutocapitalization(.never)
                    .padding(12).background(Theme.panel, in: RoundedRectangle(cornerRadius: 6)).overlay(RoundedRectangle(cornerRadius: 6).stroke(Theme.line))
                Button { searching = true; Task { let ok = await e.resolve(name); searching = false; msg = ok ? L("Found.") : L("Not found in SIMBAD.") } } label: {
                    if searching { ProgressView().frame(width: 44) } else { Image(systemName: "magnifyingglass").frame(width: 44, height: 44) }
                }.background(Theme.amberDim, in: RoundedRectangle(cornerRadius: 6)).foregroundStyle(Theme.amber)
            }
            if !msg.isEmpty { Tag(msg, color: msg == L("Found.") ? Theme.green : Theme.red) }
            HStack(spacing: 10) {
                field(L("RA  h m s"), $e.ra); field("DEC  ° ′ ″", $e.dec)
            }
            Spacer()
            Button(L("Continue")) { e.needCenter = false; e.steps = Engine.freshSteps(); e.start { await e.runPipeline() } }
                .buttonStyle(PrimaryButtonStyle()).disabled(e.ra.isEmpty || e.dec.isEmpty).opacity(e.ra.isEmpty || e.dec.isEmpty ? 0.4 : 1)
        }
        .padding(24).background(Theme.ink).presentationDetents([.medium, .large]).preferredColorScheme(.dark)
    }
    func field(_ label: String, _ b: Binding<String>) -> some View {
        VStack(alignment: .leading, spacing: 6) { Tag(label); TextField("—", text: b).font(Theme.mono(16)).padding(12)
            .background(Theme.panel, in: RoundedRectangle(cornerRadius: 6)).overlay(RoundedRectangle(cornerRadius: 6).stroke(Theme.line)) }
    }
}

// MARK: - settings and equipment profiles
struct SettingsSheet: View {
    @EnvironmentObject var e: Engine
    @Environment(\.dismiss) var dismiss
    @State private var newName = ""
    @AppStorage("lang") private var lang = "system"
    var body: some View {
        NavigationStack {
            Form {
                Section {
                    ForEach(e.profiles, id: \.name) { p in
                        Button { e.apply(profile: p) } label: {
                            HStack {
                                VStack(alignment: .leading, spacing: 2) {
                                    Text(p.name).font(Theme.mono(14, .semibold)).foregroundStyle(Theme.text)
                                    Text(p.summary).font(Theme.mono(11)).foregroundStyle(Theme.dim)
                                }
                                Spacer()
                                if p.name == e.profileName { Image(systemName: "checkmark").foregroundStyle(Theme.amber) }
                            }
                        }
                    }.onDelete { e.deleteProfiles(at: $0) }
                    HStack {
                        TextField(L("Name (e.g. C11 + ASI533)"), text: $newName).font(Theme.mono(13))
                        Button(L("Save current")) { e.saveProfile(named: newName); newName = "" }.disabled(newName.isEmpty)
                    }
                } header: { Tag(L("Equipment")) } footer: { Text(L("Profiles store focal length, pixel size, mirroring and MPC code.")).font(.caption) }
                Section {
                    Picker(L("Language"), selection: $lang) {
                        Text(L("System")).tag("system"); Text("English").tag("en"); Text("Português").tag("pt")
                    }
                } header: { Tag(L("Language")) }
                SchemaSettingsView()
            }
            .scrollContentBackground(.hidden).background(Theme.ink)
            .navigationTitle(L("Settings")).navigationBarTitleDisplayMode(.inline)
            .toolbar { ToolbarItem(placement: .confirmationAction) { Button("OK") { dismiss() } } }
        }.preferredColorScheme(.dark)
    }
}

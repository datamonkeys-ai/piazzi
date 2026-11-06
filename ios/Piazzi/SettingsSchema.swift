import SwiftUI

/// Settings schema from the core (astro_settings_schema): the settings screen is generated from it.
struct SettingsSchema: Decodable {
    struct Group: Decodable, Identifiable { let id: String, en: String, pt: String; let collapsed: Bool
        var label: String { Lang.isPT ? pt : en } }
    struct Option: Decodable, Hashable { let value: String, en: String, pt: String
        var label: String { Lang.isPT ? pt : en } }
    struct Field: Decodable, Identifiable {
        let key: String, group: String, type: String, en: String, pt: String, hint_en: String, hint_pt: String, `default`: String
        let min: Double?, max: Double?, options: [Option]?
        var id: String { key }
        var label: String { Lang.isPT ? pt : en }
        var hint: String { Lang.isPT ? hint_pt : hint_en }
    }
    let groups: [Group], settings: [Field]

    static let shared: SettingsSchema = {
        let json = String(cString: astro_settings_schema())
        return (try? JSONDecoder().decode(SettingsSchema.self, from: Data(json.utf8))) ?? SettingsSchema(groups: [], settings: [])
    }()
}

extension Engine {
    /// Resets all schema settings to the core defaults and reapplies the current equipment profile.
    func restoreDefaults() {
        for f in SettingsSchema.shared.settings { set(f.key, f.default) }
        let saved = UserDefaults.standard.string(forKey: "profile") ?? "Default"
        profileName = saved
        if let p = profiles.first(where: { $0.name == saved }) { apply(profile: p) }
    }
}

/// Form generated from the schema: one collapsible group per category ("Advanced" collapsed by default).
struct SchemaSettingsView: View {
    @EnvironmentObject var e: Engine
    @State private var open: [String: Bool] = [:]
    let schema = SettingsSchema.shared
    var body: some View {
        ForEach(schema.groups) { g in
            Section {
                DisclosureGroup(isExpanded: Binding(get: { open[g.id] ?? !g.collapsed }, set: { open[g.id] = $0 })) {
                    ForEach(schema.settings.filter { $0.group == g.id }) { f in FieldRow(f: f) }
                } label: { Tag(g.label) }
            }
        }
        Section {
            Button(L("Restore defaults")) { e.restoreDefaults() }.foregroundStyle(Theme.amber)
        }
    }
}

struct FieldRow: View {
    @EnvironmentObject var e: Engine
    let f: SettingsSchema.Field
    @State private var text = ""
    @State private var error: String?
    var value: String { e.get(f.key) }
    var body: some View {
        VStack(alignment: .leading, spacing: 4) {
            HStack {
                Text(f.label).font(.system(size: 14))
                Spacer()
                control
            }
            Text(f.hint).font(.caption).foregroundStyle(Theme.dim)
            if let err = error { Text(err).font(.caption).foregroundStyle(Theme.red) }
            if value != f.default && !(value.isEmpty && f.default.isEmpty) {
                Button(L("Restore default") + " (\(display(f.default)))") { e.set(f.key, f.default); text = f.default; error = nil }
                    .font(.caption).foregroundStyle(Theme.amber).buttonStyle(.borderless)
            }
        }
        .onAppear { text = value }
        .onChange(of: value) { _, v in if error == nil { text = v } }
    }

    func display(_ v: String) -> String {
        if let o = f.options?.first(where: { $0.value == v }) { return o.label }
        if f.type == "bool" { return v == "0" || v.isEmpty ? "✕" : "✓" }
        return v.isEmpty ? "—" : v
    }

    @ViewBuilder var control: some View {
        switch f.type {
        case "bool":
            Toggle("", isOn: Binding(get: { (Double(value) ?? 0) != 0 }, set: { e.set(f.key, $0 ? "1" : "0") })).labelsHidden()
        case "enum":
            Picker("", selection: Binding(get: { matchOption(value) }, set: { e.set(f.key, $0) })) {
                ForEach(f.options ?? [], id: \.value) { o in Text(o.label).tag(o.value) }
            }.labelsHidden().pickerStyle(.menu)
        case "number", "int":
            TextField(f.default, text: $text)
                .keyboardType(.numbersAndPunctuation).multilineTextAlignment(.trailing).font(Theme.mono(14))
                .foregroundStyle(error == nil ? Theme.amber : Theme.red).frame(maxWidth: 120)
                .onChange(of: text) { _, t in validate(t) }
        default:
            TextField(f.default, text: Binding(get: { value }, set: { e.set(f.key, $0) }))
                .multilineTextAlignment(.trailing).font(Theme.mono(14)).foregroundStyle(Theme.amber).frame(maxWidth: 160)
                .textInputAutocapitalization(.never)
        }
    }

    /// numeric enum: "7" and "7.0" are the same option
    func matchOption(_ v: String) -> String {
        let opts = f.options ?? []
        if opts.contains(where: { $0.value == v }) { return v }
        if let d = Double(v), let o = opts.first(where: { Double($0.value) == d }) { return o.value }
        return opts.first(where: { $0.value == f.default })?.value ?? v
    }

    /// Applies only valid numeric values within [min, max]; otherwise shows the error and keeps the previous value.
    func validate(_ t: String) {
        let s = t.trimmingCharacters(in: .whitespaces).replacingOccurrences(of: ",", with: ".")
        if s == value { error = nil; return }
        guard let d = Double(s), d.isFinite, f.type != "int" || d == d.rounded() else { error = L("Not a number"); return }
        if let lo = f.min, let hi = f.max, d < lo || d > hi {
            error = L("Must be between %@ and %@", fmt(lo) as NSString, fmt(hi) as NSString); return
        }
        error = nil; e.set(f.key, s)
    }
    func fmt(_ x: Double) -> String { x == x.rounded() && abs(x) < 1e9 ? String(Int(x)) : String(x) }
}

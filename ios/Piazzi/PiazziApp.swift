import SwiftUI

@main
struct PiazziApp: App {
    @StateObject private var engine = Engine()
    var body: some Scene {
        WindowGroup { ContentView().environmentObject(engine) }
    }
}

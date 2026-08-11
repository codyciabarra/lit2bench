// Lit2Bench.swift -- the native macOS app.
//
// Lit2Bench's UI is a Shiny app, so a local R server is unavoidable. What *is*
// avoidable is making the user look at a browser: this is a normal Cocoa app
// with its own NSWindow, its own Dock icon and its own menu bar, and the UI is
// hosted in a WKWebView inside it. No tabs, no address bar, no "which window was
// my app again", and ⌘Q actually quits the whole thing.
//
// Sequence on launch:
//   1. Build the menu bar and window.
//   2. Delete any stale setup page, then start launcher/bootstrap.sh as a child
//      process with LIT2BENCH_NATIVE=1 (so it doesn't open a browser).
//   3. Poll for the setup page bootstrap.sh writes, and load it.
//   4. Do nothing further. The setup page already reloads itself while work is
//      happening and redirects to the server when it's up, so the WKWebView
//      follows it into the app on its own -- the same mechanism that used to
//      drive a browser tab, now driving a native window.
//
// Hosting a web UI natively means re-supplying the things a browser gave us for
// free, and each of these is load-bearing for a specific tool:
//   * an Edit menu with the standard selectors, or ⌘C/⌘V don't work at all;
//   * runOpenPanelWith, or <input type="file"> does nothing -- that's the
//     Cryptic Engine's BAM uploads and Plasmid QC's read files;
//   * WKDownload, or every export button in the app silently fails;
//   * link interception, so github.com opens in Safari rather than replacing
//     the app UI with a web page the user can't navigate back from.
//
// Built by installer/macos/build.sh as a universal binary. Single file on
// purpose: swiftc compiles it directly, with no Xcode project to maintain.

import Cocoa
import WebKit

// MARK: - Paths

/// Everything writable lives beside the notebook, under Application Support --
/// the bundle itself is read-only. Mirrors R/paths.R and bootstrap.sh.
let supportDir: URL = {
    let base = FileManager.default.urls(for: .applicationSupportDirectory, in: .userDomainMask)[0]
    return base.appendingPathComponent("Lit2Bench", isDirectory: true)
}()

let statusPageURL = supportDir
    .appendingPathComponent("setup", isDirectory: true)
    .appendingPathComponent("status.html")

let bootstrapURL = Bundle.main.resourceURL!
    .appendingPathComponent("launcher", isDirectory: true)
    .appendingPathComponent("bootstrap.sh")

// MARK: - App

final class AppDelegate: NSObject, NSApplicationDelegate {
    private var window: NSWindow!
    private var web: WKWebView!
    private var server: Process?
    private var pollTimer: Timer?
    /// WKDownload doesn't hand the destination back on completion, so keep it.
    private var downloadDestinations: [ObjectIdentifier: URL] = [:]

    func applicationDidFinishLaunching(_ notification: Notification) {
        buildMenu()
        buildWindow()
        startServer()
        waitForSetupPage()
        NSApp.activate(ignoringOtherApps: true)
    }

    func applicationShouldTerminateAfterLastWindowClosed(_ app: NSApplication) -> Bool { true }

    func applicationWillTerminate(_ notification: Notification) {
        pollTimer?.invalidate()
        // SIGTERM, not SIGKILL: bootstrap.sh traps it and takes the R server
        // down with it. Killing outright would orphan R still holding the port.
        server?.terminate()
        server?.waitUntilExit()
    }

    // MARK: Window

    private func buildWindow() {
        let config = WKWebViewConfiguration()
        // The R server is same-origin for everything the UI does, but the setup
        // page is a file:// document that navigates to http://127.0.0.1 -- allow it.
        config.preferences.setValue(true, forKey: "allowFileAccessFromFileURLs")
        config.defaultWebpagePreferences.allowsContentJavaScript = true

        web = WKWebView(frame: .zero, configuration: config)
        web.navigationDelegate = self
        web.uiDelegate = self
        web.allowsBackForwardNavigationGestures = false
        // Adapts to light/dark, so neither theme flashes the wrong colour while
        // the first paint is still pending.
        web.underPageBackgroundColor = .windowBackgroundColor
        if web.responds(to: Selector(("setInspectable:"))) {
            web.setValue(true, forKey: "inspectable")   // Web Inspector, for debugging a tool
        }

        window = NSWindow(
            contentRect: NSRect(x: 0, y: 0, width: 1440, height: 900),
            styleMask: [.titled, .closable, .miniaturizable, .resizable, .fullSizeContentView],
            backing: .buffered, defer: false)
        window.title = "Lit2Bench"
        window.titlebarAppearsTransparent = false
        window.minSize = NSSize(width: 1000, height: 680)
        window.contentView = web
        window.setFrameAutosaveName("Lit2BenchMainWindow")
        window.center()
        window.makeKeyAndOrderFront(nil)
    }

    // MARK: Server

    private func startServer() {
        // A stale page from a previous run would say "ready" and redirect to a
        // port nothing is listening on any more.
        try? FileManager.default.removeItem(at: statusPageURL)

        let p = Process()
        p.executableURL = URL(fileURLWithPath: "/bin/bash")
        p.arguments = [bootstrapURL.path]
        var env = ProcessInfo.processInfo.environment
        env["LIT2BENCH_NATIVE"] = "1"     // don't open a browser; we are the window
        p.environment = env
        // bootstrap.sh logs to files; anything on the pipe would eventually
        // fill the buffer and wedge it.
        p.standardOutput = FileHandle.nullDevice
        p.standardError = FileHandle.nullDevice

        do {
            try p.run()
            server = p
        } catch {
            showFailure("Couldn't start Lit2Bench's setup process.\n\n\(error.localizedDescription)")
        }
    }

    /// bootstrap.sh writes the setup page as its first action, so this normally
    /// resolves on the first or second tick.
    private func waitForSetupPage() {
        var waited = 0.0
        pollTimer = Timer.scheduledTimer(withTimeInterval: 0.08, repeats: true) { [weak self] timer in
            guard let self else { timer.invalidate(); return }
            if FileManager.default.fileExists(atPath: statusPageURL.path) {
                timer.invalidate()
                self.web.loadFileURL(statusPageURL, allowingReadAccessTo: supportDir)
                return
            }
            waited += 0.08
            if waited > 20 {
                timer.invalidate()
                self.showFailure("Lit2Bench's setup process didn't start.\n\nCheck \(supportDir.path)/logs/setup.log")
            }
        }
    }

    private func showFailure(_ message: String) {
        let alert = NSAlert()
        alert.alertStyle = .critical
        alert.messageText = "Lit2Bench couldn't start"
        alert.informativeText = message
        alert.addButton(withTitle: "Quit")
        alert.addButton(withTitle: "Open Log Folder")
        if alert.runModal() == .alertSecondButtonReturn {
            NSWorkspace.shared.open(supportDir.appendingPathComponent("logs"))
        }
        NSApp.terminate(nil)
    }

    // MARK: Menu actions

    @objc private func reload() {
        // The setup page is a plain file with no server behind it; reloading the
        // app itself is what's meant here.
        if web.url?.isFileURL == true { web.reload() } else { web.reloadFromOrigin() }
    }
    @objc private func zoomIn()     { web.pageZoom = min(web.pageZoom + 0.1, 3.0) }
    @objc private func zoomOut()    { web.pageZoom = max(web.pageZoom - 0.1, 0.5) }
    @objc private func zoomActual() { web.pageZoom = 1.0 }

    @objc private func openInBrowser() {
        // Escape hatch: a couple of things (printing a big figure, say) are just
        // nicer in a real browser.
        if let url = web.url, !url.isFileURL { NSWorkspace.shared.open(url) }
    }

    @objc private func openLogs() {
        NSWorkspace.shared.open(supportDir.appendingPathComponent("logs"))
    }

    @objc private func openDataFolder() {
        NSWorkspace.shared.open(supportDir)
    }

    private func buildMenu() {
        let mainMenu = NSMenu()

        // Application menu. The title is ignored by macOS -- it always shows the
        // bundle name -- but the items are not.
        let appItem = NSMenuItem()
        let appMenu = NSMenu()
        appMenu.addItem(withTitle: "About Lit2Bench", action: #selector(NSApplication.orderFrontStandardAboutPanel(_:)), keyEquivalent: "")
        appMenu.addItem(.separator())
        appMenu.addItem(withTitle: "Show Data Folder", action: #selector(openDataFolder), keyEquivalent: "")
        appMenu.addItem(withTitle: "Show Logs", action: #selector(openLogs), keyEquivalent: "")
        appMenu.addItem(.separator())
        appMenu.addItem(withTitle: "Hide Lit2Bench", action: #selector(NSApplication.hide(_:)), keyEquivalent: "h")
        let hideOthers = appMenu.addItem(withTitle: "Hide Others", action: #selector(NSApplication.hideOtherApplications(_:)), keyEquivalent: "h")
        hideOthers.keyEquivalentModifierMask = [.command, .option]
        appMenu.addItem(withTitle: "Show All", action: #selector(NSApplication.unhideAllApplications(_:)), keyEquivalent: "")
        appMenu.addItem(.separator())
        appMenu.addItem(withTitle: "Quit Lit2Bench", action: #selector(NSApplication.terminate(_:)), keyEquivalent: "q")
        appItem.submenu = appMenu
        mainMenu.addItem(appItem)

        let fileItem = NSMenuItem()
        let fileMenu = NSMenu(title: "File")
        let browser = fileMenu.addItem(withTitle: "Open in Browser", action: #selector(openInBrowser), keyEquivalent: "b")
        browser.keyEquivalentModifierMask = [.command, .shift]
        fileMenu.addItem(.separator())
        fileMenu.addItem(withTitle: "Close Window", action: #selector(NSWindow.performClose(_:)), keyEquivalent: "w")
        fileItem.submenu = fileMenu
        mainMenu.addItem(fileItem)

        // Without this menu, ⌘C and ⌘V do nothing anywhere in the app: the
        // responder chain has no cut/copy/paste entry to route them through.
        let editItem = NSMenuItem()
        let editMenu = NSMenu(title: "Edit")
        editMenu.addItem(withTitle: "Undo", action: Selector(("undo:")), keyEquivalent: "z")
        let redo = editMenu.addItem(withTitle: "Redo", action: Selector(("redo:")), keyEquivalent: "z")
        redo.keyEquivalentModifierMask = [.command, .shift]
        editMenu.addItem(.separator())
        editMenu.addItem(withTitle: "Cut", action: #selector(NSText.cut(_:)), keyEquivalent: "x")
        editMenu.addItem(withTitle: "Copy", action: #selector(NSText.copy(_:)), keyEquivalent: "c")
        editMenu.addItem(withTitle: "Paste", action: #selector(NSText.paste(_:)), keyEquivalent: "v")
        editMenu.addItem(withTitle: "Select All", action: #selector(NSText.selectAll(_:)), keyEquivalent: "a")
        editItem.submenu = editMenu
        mainMenu.addItem(editItem)

        let viewItem = NSMenuItem()
        let viewMenu = NSMenu(title: "View")
        viewMenu.addItem(withTitle: "Reload", action: #selector(reload), keyEquivalent: "r")
        viewMenu.addItem(.separator())
        viewMenu.addItem(withTitle: "Actual Size", action: #selector(zoomActual), keyEquivalent: "0")
        viewMenu.addItem(withTitle: "Zoom In", action: #selector(zoomIn), keyEquivalent: "+")
        viewMenu.addItem(withTitle: "Zoom Out", action: #selector(zoomOut), keyEquivalent: "-")
        viewMenu.addItem(.separator())
        let full = viewMenu.addItem(withTitle: "Enter Full Screen", action: #selector(NSWindow.toggleFullScreen(_:)), keyEquivalent: "f")
        full.keyEquivalentModifierMask = [.command, .control]
        viewItem.submenu = viewMenu
        mainMenu.addItem(viewItem)

        let windowItem = NSMenuItem()
        let windowMenu = NSMenu(title: "Window")
        windowMenu.addItem(withTitle: "Minimize", action: #selector(NSWindow.performMiniaturize(_:)), keyEquivalent: "m")
        windowMenu.addItem(withTitle: "Zoom", action: #selector(NSWindow.performZoom(_:)), keyEquivalent: "")
        windowMenu.addItem(.separator())
        windowMenu.addItem(withTitle: "Bring All to Front", action: #selector(NSApplication.arrangeInFront(_:)), keyEquivalent: "")
        windowItem.submenu = windowMenu
        mainMenu.addItem(windowItem)

        NSApp.mainMenu = mainMenu
        NSApp.windowsMenu = windowMenu
    }
}

// MARK: - Navigation

extension AppDelegate: WKNavigationDelegate {
    /// Keep the window on the app, and hand anything else to the real browser.
    /// Without this, clicking the GitHub link in the About tab would replace the
    /// entire UI with a web page and no way back.
    func webView(_ webView: WKWebView,
                 decidePolicyFor navigationAction: WKNavigationAction,
                 decisionHandler: @escaping (WKNavigationActionPolicy) -> Void) {
        guard let url = navigationAction.request.url else { return decisionHandler(.allow) }
        let host = url.host ?? ""
        let isLocal = url.isFileURL
            || host == "127.0.0.1" || host == "localhost"
            || ["about", "data", "blob"].contains(url.scheme ?? "")
        let isUserClick = navigationAction.navigationType == .linkActivated
            || navigationAction.targetFrame == nil

        if !isLocal && isUserClick {
            NSWorkspace.shared.open(url)
            return decisionHandler(.cancel)
        }
        decisionHandler(.allow)
    }

    /// Shiny's downloadHandler sends Content-Disposition: attachment. WKWebView
    /// would happily *render* some of those (the HTML figure exports especially),
    /// so trust the header over the MIME type.
    func webView(_ webView: WKWebView,
                 decidePolicyFor navigationResponse: WKNavigationResponse,
                 decisionHandler: @escaping (WKNavigationResponsePolicy) -> Void) {
        if let http = navigationResponse.response as? HTTPURLResponse,
           let disposition = http.value(forHTTPHeaderField: "Content-Disposition"),
           disposition.lowercased().contains("attachment") {
            return decisionHandler(.download)
        }
        decisionHandler(navigationResponse.canShowMIMEType ? .allow : .download)
    }

    func webView(_ webView: WKWebView, navigationResponse: WKNavigationResponse, didBecome download: WKDownload) {
        download.delegate = self
    }

    func webView(_ webView: WKWebView, navigationAction: WKNavigationAction, didBecome download: WKDownload) {
        download.delegate = self
    }

    func webView(_ webView: WKWebView, didFail navigation: WKNavigation!, withError error: Error) {
        reportNavigationFailure(error)
    }

    func webView(_ webView: WKWebView, didFailProvisionalNavigation navigation: WKNavigation!, withError error: Error) {
        reportNavigationFailure(error)
    }

    private func reportNavigationFailure(_ error: Error) {
        // Shiny drops the connection when the server is restarting or has been
        // stopped; a cancelled load is normal and not worth a dialog.
        let ns = error as NSError
        guard ns.domain != NSURLErrorDomain || ns.code != NSURLErrorCancelled else { return }
        window.title = "Lit2Bench — disconnected"
    }
}

// MARK: - UI delegate (file pickers, dialogs, window.open)

extension AppDelegate: WKUIDelegate {
    /// <input type="file"> is inert in a WKWebView until this exists. The BAM
    /// uploads in the Cryptic Engine and the read files in Plasmid QC both need it.
    func webView(_ webView: WKWebView,
                 runOpenPanelWith parameters: WKOpenPanelParameters,
                 initiatedByFrame frame: WKFrameInfo,
                 completionHandler: @escaping ([URL]?) -> Void) {
        let panel = NSOpenPanel()
        panel.allowsMultipleSelection = parameters.allowsMultipleSelection
        panel.canChooseDirectories = parameters.allowsDirectories
        panel.canChooseFiles = true
        panel.resolvesAliases = true
        panel.beginSheetModal(for: window) { response in
            completionHandler(response == .OK ? panel.urls : nil)
        }
    }

    func webView(_ webView: WKWebView, createWebViewWith configuration: WKWebViewConfiguration,
                 for navigationAction: WKNavigationAction,
                 windowFeatures: WKWindowFeatures) -> WKWebView? {
        if let url = navigationAction.request.url { NSWorkspace.shared.open(url) }
        return nil   // never spawn a second, chromeless window
    }

    func webView(_ webView: WKWebView, runJavaScriptAlertPanelWithMessage message: String,
                 initiatedByFrame frame: WKFrameInfo, completionHandler: @escaping () -> Void) {
        let alert = NSAlert()
        alert.messageText = "Lit2Bench"
        alert.informativeText = message
        alert.beginSheetModal(for: window) { _ in completionHandler() }
    }

    func webView(_ webView: WKWebView, runJavaScriptConfirmPanelWithMessage message: String,
                 initiatedByFrame frame: WKFrameInfo, completionHandler: @escaping (Bool) -> Void) {
        let alert = NSAlert()
        alert.messageText = "Lit2Bench"
        alert.informativeText = message
        alert.addButton(withTitle: "OK")
        alert.addButton(withTitle: "Cancel")
        alert.beginSheetModal(for: window) { completionHandler($0 == .alertFirstButtonReturn) }
    }
}

// MARK: - Downloads

extension AppDelegate: WKDownloadDelegate {
    func download(_ download: WKDownload,
                  decideDestinationUsing response: URLResponse,
                  suggestedFilename: String,
                  completionHandler: @escaping (URL?) -> Void) {
        let downloads = FileManager.default.urls(for: .downloadsDirectory, in: .userDomainMask)[0]
        let name = suggestedFilename.isEmpty ? "lit2bench-export" : suggestedFilename

        // WKDownload refuses to overwrite, so pick the first free "name (n).ext"
        // rather than failing the export on a repeat run.
        var destination = downloads.appendingPathComponent(name)
        if FileManager.default.fileExists(atPath: destination.path) {
            let ext = destination.pathExtension
            let stem = destination.deletingPathExtension().lastPathComponent
            var n = 2
            repeat {
                let candidate = ext.isEmpty ? "\(stem) (\(n))" : "\(stem) (\(n)).\(ext)"
                destination = downloads.appendingPathComponent(candidate)
                n += 1
            } while FileManager.default.fileExists(atPath: destination.path) && n < 1000
        }

        downloadDestinations[ObjectIdentifier(download)] = destination
        completionHandler(destination)
    }

    func downloadDidFinish(_ download: WKDownload) {
        let key = ObjectIdentifier(download)
        if let url = downloadDestinations.removeValue(forKey: key) {
            // Same feedback a browser gives: bounce the Dock and reveal the file.
            NSWorkspace.shared.activateFileViewerSelecting([url])
        }
    }

    func download(_ download: WKDownload, didFailWithError error: Error, resumeData: Data?) {
        downloadDestinations.removeValue(forKey: ObjectIdentifier(download))
        let alert = NSAlert()
        alert.alertStyle = .warning
        alert.messageText = "Download failed"
        alert.informativeText = error.localizedDescription
        alert.beginSheetModal(for: window, completionHandler: nil)
    }
}

// MARK: - Entry point

let app = NSApplication.shared
let delegate = AppDelegate()
app.delegate = delegate
app.setActivationPolicy(.regular)
app.run()

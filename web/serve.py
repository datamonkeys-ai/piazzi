#!/usr/bin/env python3
"""Local server for Piazzi web: static files + proxy to the astronomy services.

The browser cannot call VizieR/SkyBoT/SIMBAD directly (CORS); this server proxies only to allowed hosts.

Usage: python3 web/serve.py [port]
Environment variables (see .env.example): PORT (8777), HOST (127.0.0.1; 0.0.0.0 exposes on the LAN), PROXY_TIMEOUT (90 s), LOG_LEVEL (INFO),
EXTRA_HOSTS (extra names/IPs accepted in the Host header, comma-separated — e.g. the Docker host IP).
The Host header is validated (localhost/127.0.0.1/[::1] and, when not bound to 127.0.0.1, this machine's IPs) against DNS rebinding.
Routes: /            static files from web/
        /proxy?u=URL proxy to a host in the ALLOWED list
        /healthz     health check (used by the Docker healthcheck)
"""

import http.server
import json
import logging
import os
import socket
import sys
import urllib.error
import urllib.parse
import urllib.request

ALLOWED = frozenset(
    {
        "vizier.idia.ac.za",
        "vizier.cds.unistra.fr",
        "vizier.cfa.harvard.edu",
        "ssp.imcce.fr",
        "simbad.cfa.harvard.edu",
        "simbad.cds.unistra.fr",
        "cds.unistra.fr",
    }
)
ROOT = os.path.dirname(os.path.abspath(__file__))
PROXY_TIMEOUT = float(os.environ.get("PROXY_TIMEOUT", "90"))
MAX_BYTES = int(os.environ.get("PROXY_MAX_BYTES", str(64 * 1024 * 1024)))  # larger responses are rejected
log = logging.getLogger("piazzi.web")
# Security headers on every response (no third-party framing, no content-type sniffing).
CABECALHOS_SEG = {
    "X-Content-Type-Options": "nosniff",
    "X-Frame-Options": "DENY",
    "Content-Security-Policy": "frame-ancestors 'none'; object-src 'none'; base-uri 'none'",
    "Referrer-Policy": "no-referrer",
}
HOSTS_OK: set[str] = {"localhost", "127.0.0.1", "[::1]"}  # names accepted in Host (without port); see configurar_hosts()


def ips_proprios() -> set[str]:
    ips = set()
    try:
        ips.update(socket.gethostbyname_ex(socket.gethostname())[2])
    except OSError:
        pass
    ip = lan_ip()
    if ip != "localhost":
        ips.add(ip)
    return ips


def configurar_hosts(host: str) -> None:
    """Outside local mode also accepts this machine's IPs (and EXTRA_HOSTS); never arbitrary names."""
    if host != "127.0.0.1":
        HOSTS_OK.update(ips_proprios())
        if host not in ("", "0.0.0.0", "::"):
            HOSTS_OK.add(host)
    HOSTS_OK.update(h.strip().lower() for h in os.environ.get("EXTRA_HOSTS", "").split(",") if h.strip())


def host_ok(cab: str) -> bool:
    """Compares only the name in the Host header (the port changes with Docker port mapping)."""
    cab = cab.strip().lower()
    nome = cab.rsplit(":", 1)[0] if not cab.endswith("]") else cab
    return nome in HOSTS_OK


def allowed_url(url: str) -> bool:
    """HTTPS only and only to listed hosts; prevents use as an open proxy."""
    p = urllib.parse.urlparse(url)
    return p.scheme == "https" and (p.hostname or "") in ALLOWED


class SafeRedirect(urllib.request.HTTPRedirectHandler):
    """Revalidates the target of every redirect (otherwise an allowed host could lead to any other)."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if not allowed_url(newurl):
            raise urllib.error.HTTPError(newurl, 403, "redirect to a disallowed host", headers, fp)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


OPENER = urllib.request.build_opener(SafeRedirect)


class Handler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *a, **k):
        super().__init__(*a, directory=ROOT, **k)

    def end_headers(self):
        for k, v in CABECALHOS_SEG.items():
            self.send_header(k, v)
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def _send(self, code: int, body: bytes, ctype: str) -> None:
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_HEAD(self):
        if not host_ok(self.headers.get("Host", "")):
            return self._send(403, b"", "text/plain; charset=utf-8")
        return super().do_HEAD()

    def do_GET(self):
        if self.path == "/healthz":
            ok = os.path.exists(os.path.join(ROOT, "astro.wasm"))
            body = json.dumps({"status": "ok" if ok else "degraded", "wasm": ok}).encode()
            return self._send(200 if ok else 503, body, "application/json")
        if not host_ok(self.headers.get("Host", "")):  # DNS rebinding: unknown name in Host
            log.warning("Host rejected: %s", self.headers.get("Host", "")[:80])
            return self._send(403, b"host rejected", "text/plain; charset=utf-8")
        if not self.path.startswith("/proxy?"):
            return super().do_GET()
        url = urllib.parse.parse_qs(urllib.parse.urlparse(self.path).query).get("u", [""])[0]
        if not allowed_url(url):
            log.warning("proxy rejected: %s", url[:140])
            self.send_error(403, "host not allowed")
            return
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Piazzi-web/1.0"})
            with OPENER.open(req, timeout=PROXY_TIMEOUT) as r:
                data = r.read(MAX_BYTES + 1)
                if len(data) > MAX_BYTES:
                    log.warning("proxy %s: response exceeds %d bytes", url[:140], MAX_BYTES)
                    self.send_error(502, "response too large")
                    return
                self._send(r.status, data, r.headers.get("Content-Type", "application/octet-stream"))
            log.info("proxy %s -> %s", urllib.parse.urlparse(url).hostname, r.status)
        except urllib.error.HTTPError as e:
            log.warning("proxy %s -> HTTP %s", url[:140], e.code)
            self.send_error(e.code, str(e.reason))
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            log.error("proxy %s failed: %s", url[:140], e)
            self.send_error(502, str(e)[:200])

    def log_message(self, fmt, *args):  # static file access is logged at DEBUG; the proxy already logs above
        log.debug(fmt, *args)


def lan_ip() -> str:
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as s:
            s.connect(("192.0.2.1", 80))  # documentation address: sends no packets, only selects the interface
            return s.getsockname()[0]
    except OSError:
        return "localhost"


def main() -> None:
    logging.basicConfig(level=os.environ.get("LOG_LEVEL", "INFO"), format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    port = int(sys.argv[1]) if len(sys.argv) > 1 else int(os.environ.get("PORT", "8777"))
    host = os.environ.get("HOST", "127.0.0.1")
    configurar_hosts(host)
    log.info("Piazzi web: http://localhost:%d  ·  on the LAN (iPhone): http://%s:%d", port, lan_ip(), port)
    http.server.ThreadingHTTPServer((host, port), Handler).serve_forever()


if __name__ == "__main__":
    main()

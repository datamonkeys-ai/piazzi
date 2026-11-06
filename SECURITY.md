# Security

## Reporting a vulnerability
Use GitHub private vulnerability reporting (Security tab → Report a vulnerability) with steps to reproduce. Do not open a public issue before the fix.
Initial response within 7 days.

## Threat model (summary)
| Surface | Risk | Mitigation |
|---|---|---|
| FITS/SBIG files opened by the user | out-of-bounds read, malformed file | readers validate dimensions, BITPIX and size before allocating; errors become messages, not crashes |
| External processes (curl, open/xdg-open) in the desktop app | command injection via e-mail or path | execution without a shell (`posix_spawnp` / `CreateProcess` with escaped arguments); e-mail validated |
| `web/serve.py` (CORS proxy) | open proxy / SSRF / DNS rebinding | `https` only, to a fixed list of astronomy hosts; redirects revalidated; Host validated; listens on 127.0.0.1 by default |
| Submission to the MPC | wrong report in the observer's name | explicit confirmation; test mode checked by default |
| Dependencies and CI | vulnerable or tampered components | dependency-free core; actions and Docker images pinned by SHA/digest; Dependabot |

The app requires and stores no credentials.

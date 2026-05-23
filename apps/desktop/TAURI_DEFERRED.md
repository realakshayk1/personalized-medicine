Tauri shell deferred per Wave 1 decision (Windows toolchain risk — Rust/MSVC setup adds friction to the Day-14 demo path).
The browser-served Next.js app is the Day-14 demo target; all UI and logic runs in the browser via `next dev` or `next start`.
A future session will wrap this Next.js app in a Tauri v2 shell once the Windows Rust toolchain is verified.
No `src-tauri/` directory exists in this workspace; do not create it until that session begins.

// Piazzi web — same C++ core (WebAssembly) as the iPhone app and the CLI.
'use strict';
const $ = (s) => document.querySelector(s);
const el = (t, cls, html) => { const e = document.createElement(t); if (cls) e.className = cls; if (html !== undefined) e.innerHTML = html; return e; };
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const C = { amber: '#ffb547', red: '#ff5c4f', cyan: '#80d4ff', green: '#9be564', white: '#ffffff' };
const proxy = (u) => 'proxy?u=' + encodeURIComponent(u); // relative: works at / (serve.py) and under /piazzi/app/

// ---------------------------------------------------------------- background: star chart
function drawChart(cv, seed = 7, vignette = true) {
  const dpr = devicePixelRatio || 1, w = cv.clientWidth || innerWidth, h = cv.clientHeight || innerHeight;
  cv.width = w * dpr; cv.height = h * dpr; const g = cv.getContext('2d'); g.scale(dpr, dpr); g.clearRect(0, 0, w, h);
  g.strokeStyle = 'rgba(255,255,255,.035)'; g.lineWidth = .6;
  for (let i = 0; i <= 8; i++) { const x = w * i / 8; g.beginPath(); g.moveTo(x, 0); g.quadraticCurveTo(x + (x - w / 2) * .03, h / 2, x + (x - w / 2) * .12, h); g.stroke(); }
  for (let j = 0; j <= 12; j++) { const y = h * j / 12; g.beginPath(); g.moveTo(0, y); g.quadraticCurveTo(w / 2, y - 18, w, y); g.stroke(); }
  let s = BigInt(seed); const rnd = () => { s = (s * 6364136223846793005n + 1442695040888963407n) & 0xffffffffffffffffn; return Number(s >> 33n) / 2147483648; };
  for (let i = 0; i < 260; i++) { const x = rnd() * w, y = rnd() * h, m = rnd(), r = .35 + Math.pow(m, 6) * 2.4;
    g.fillStyle = `rgba(255,255,255,${.15 + m * .5})`; g.beginPath(); g.arc(x, y, r / 2, 0, 7); g.fill(); }
  if (vignette) { const gr = g.createRadialGradient(w / 2, h * .35, Math.min(w, h) * .25, w / 2, h * .35, Math.max(w, h) * .8);
    gr.addColorStop(0, 'rgba(0,0,0,0)'); gr.addColorStop(1, 'rgba(0,0,0,.65)'); g.fillStyle = gr; g.fillRect(0, 0, w, h); }
}
function reticle(g, x, y, r, gap = .45) {
  const l = r * gap; g.beginPath();
  for (const [sx, sy] of [[-1, -1], [1, -1], [-1, 1], [1, 1]]) { const cx = x + r * sx, cy = y + r * sy; g.moveTo(cx, cy - l * sy); g.lineTo(cx, cy); g.lineTo(cx - l * sx, cy); }
  g.stroke();
}

// ---------------------------------------------------------------- core (WASM)
const A = { M: null, s: 0, images: [], catalog: null, ra: '', dec: '', object: '', result: null, known: [], identified: {}, steps: [], status: '',
  profile: lsGet('profile') || '', task: 0, example: false, schema: null };
const PROFILES_BUILTIN = [];
function lsGet(k) { try { return localStorage.getItem(k); } catch (e) { return null; } }
function lsSet(k, v) { try { localStorage.setItem(k, v); } catch (e) { /* no storage */ } }
const profiles = () => { try { return PROFILES_BUILTIN.concat(JSON.parse(lsGet('profiles') || '[]')); } catch (e) { return PROFILES_BUILTIN; } };
const profileName = () => A.profile || t('default_profile');
let api;
async function boot() {
  A.M = await AstroModule();
  const w = (n, r, a) => A.M.cwrap(n, r, a);
  api = { new: w('astro_new', 'number', []), loadCfg: w('astro_load_cfg', null, ['number', 'string']), set: w('astro_set', null, ['number', 'string', 'string']),
    get: w('astro_get', 'string', ['number', 'string']), add: w('astro_add_image', 'number', ['number', 'string']), clear: w('astro_clear_images', null, ['number']),
    url: w('astro_catalog_url', 'string', ['number', 'string', 'string']), reduce: w('astro_reduce', 'number', ['number', 'string', 'string', 'string', 'string']),
    json: w('astro_result_json', 'string', ['number']), report: w('astro_report', 'string', ['number', 'number']), size: w('astro_frame_size', 'number', ['number', 'number', 'number', 'number']),
    render: w('astro_render', null, ['number', 'number', 'number']), s2p: w('astro_sky_to_pix', 'number', ['number', 'number', 'number', 'number', 'number', 'number']),
    p2s: w('astro_pix_to_sky', 'number', ['number', 'number', 'number', 'number', 'number', 'number']), center: w('astro_frame_center', 'number', ['number', 'number', 'number', 'number']),
    desig: w('astro_set_mover_desig', null, ['number', 'number', 'string']), schema: w('astro_settings_schema', 'string', []) };
  A.s = api.new(); A.M.FS.mkdir('/in');
  try { A.schema = JSON.parse(api.schema()); } catch (e) { A.schema = { groups: [], settings: [] }; }
  $('#version').textContent = `Piazzi ${A.M.ccall('astro_version', 'string', [], [])} · open source · MIT`;
  const p = profiles().find((x) => x.name === A.profile); if (p) applyProfile(p);
  $('#btn-settings').textContent = profileName() + ' ⚙︎';
  if (/[?&]example=1/.test(location.search)) startExample();
}
const set = (k, v) => api.set(A.s, k, String(v));
const get = (k) => api.get(A.s, k);
function applyProfile(p) { set('focal_length_mm', p.focal); set('pixel_width_mm', p.pw); set('pixel_height_mm', p.ph); set('flip_horizontal', p.flip); set('mpc_code', p.mpc);
  A.profile = p.name; lsSet('profile', p.name); const b = $('#btn-settings'); if (b) b.textContent = p.name + ' ⚙︎'; }
function addFile(name, bytes) {
  const path = '/in/' + name; A.M.FS.writeFile(path, bytes);
  const ext = name.split('.').pop().toLowerCase();
  if (ext === 'cfg') api.loadCfg(A.s, path); else if (ext === 'txt' || ext === 'tsv') A.catalog = path; else { A.images.push({ name, path, bytes }); api.add(A.s, path); }
}
function clearSession() { api.clear(A.s); A.images = []; A.catalog = null; A.ra = A.dec = A.object = ''; A.result = null; A.known = []; A.identified = {}; A.inspect = null; A.example = false; }
const frames = () => (A.result && A.result.frames) || [];
const movers = () => (A.result && A.result.movers) || [];
// display label: without a designation (unidentified object), the track's temporary designation (trkSub) or "object N"
const moverLabel = (m) => m.desig || m.trksub || t('object_n', movers().indexOf(m) + 1);
// Round to tenths of a second before splitting (never "59 60.0"); RA (unsigned) modulo 24.
function fmt(v, sign) {
  if (!Number.isFinite(v)) v = 0;
  let t = Math.round(Math.abs(v) * 36000), neg = v < 0;
  if (!sign) { const day = 24 * 36000; t = neg ? (day - t % day) % day : t % day; neg = false; }
  if (t === 0) neg = false;
  const p2 = (n) => String(n).padStart(2, '0'), ds = t % 600;
  return (sign ? (neg ? '-' : '+') : '') + p2(Math.floor(t / 36000)) + ' ' + p2(Math.floor(t / 600) % 60) + ' ' + p2(Math.floor(ds / 10)) + '.' + (ds % 10);
}
function sexa(s) { const p = s.trim().split(/\s+/).map(Number); if (p.length < 3) return Number(s); const sg = s.trim().startsWith('-') ? -1 : 1; return sg * (Math.abs(p[0]) + p[1] / 60 + p[2] / 3600); }
function dbl2(fn, ...args) { const p = A.M._malloc(16); try { const ok = fn(A.s, ...args, p, p + 8); return ok ? [A.M.HEAPF64[p >> 3], A.M.HEAPF64[(p >> 3) + 1]] : null; } finally { A.M._free(p); } }

// ---------------------------------------------------------------- network
async function fetchRetry(url, { tries = 3, min = 0, text = true, timeout = 60000 } = {}) {
  let last;
  for (let t = 0; t < tries; t++) {
    const ac = new AbortController(); const id = setTimeout(() => ac.abort(), timeout);
    try { const r = await fetch(proxy(url), { signal: ac.signal }); clearTimeout(id);
      if (r.ok) { const b = text ? await r.text() : new Uint8Array(await r.arrayBuffer()); if ((b.length || 0) >= min) return b; }
      last = new Error('HTTP ' + r.status);
    } catch (e) { clearTimeout(id); last = e; }
    if (t + 1 < tries) await new Promise((r) => setTimeout(r, 2000 * (t + 1)));
  }
  throw last;
}
const VIZIER = ['vizier.idia.ac.za', 'vizier.cds.unistra.fr', 'vizier.cfa.harvard.edu'];
async function skybot(jd, ra, dec, rd, obs) {
  const o = !obs || obs === 'XXX' ? '500' : obs;
  const t = await fetchRetry(`https://ssp.imcce.fr/webservices/skybot/api/conesearch.php?-ep=${jd}&-ra=${ra}&-dec=${dec}&-rd=${rd}&-mime=text&-output=object&-observer=${o}&-objFilter=110&-from=Piazzi`, { tries: 2, timeout: 25000 });
  return t.split('\n').filter((l) => !l.startsWith('#') && l.includes('|')).map((l) => { const c = l.split('|').map((x) => x.trim());
    return { number: c[0] && c[0] !== '-' ? c[0] : null, name: c[1], ra: sexa(c[2]) * 15, dec: sexa(c[3]), cls: c[4], v: Number(c[5]), sep: Number(c[7]) }; })
    .filter((a) => isFinite(a.v)).sort((a, b) => a.sep - b.sep);
}
async function vizier(source, ra, dec, r, cols) {
  for (const h of VIZIER) { try {
    const t = await fetchRetry(`https://${h}/viz-bin/asu-tsv?-source=${source}&-c=${ra}%20${dec >= 0 ? '%2B' : ''}${dec}&-c.rs=${r}&-out=${cols}&-oc.form=d&-sort=_r`, { tries: 1, timeout: 15000 });
    return t.split('\n').filter((l) => l && !l.startsWith('#')).slice(3).map((l) => l.split('\t').map((x) => x.trim()));
  } catch (e) { /* next mirror */ } }
  return null;
}
async function simbad(ra, dec, r) {
  const q = `SELECT TOP 3 main_id, otype, DISTANCE(POINT('ICRS',ra,dec),POINT('ICRS',${ra},${dec}))*3600 AS d FROM basic WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',${ra},${dec},${r / 3600}))=1 ORDER BY d`;
  for (const h of ['simbad.cfa.harvard.edu', 'simbad.cds.unistra.fr']) { try {
    const t = await fetchRetry(`https://${h}/simbad/sim-tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=${encodeURIComponent(q)}`, { tries: 1, timeout: 15000 });
    return t.split('\n').slice(1).filter(Boolean).map((l) => { const c = l.split(',').map((x) => x.replace(/"/g, '').trim()); return { id: c[0], type: c[1], sep: Number(c[2]) }; });
  } catch (e) { } }
  return null;
}

// ---------------------------------------------------------------- steps
const STEPS = ['s_images', 's_catalog', 's_plate', 's_motion', 's_ident'];
function freshSteps() { A.steps = STEPS.map((k) => ({ k, st: 'pending', d: '', at: '' })); A.t0 = Date.now(); renderSession(); }
function step(i, st, d = '') { const s = A.steps[i]; s.st = st; s.d = d; const e = Math.round((Date.now() - A.t0) / 1000); s.at = String(Math.floor(e / 60)).padStart(2, '0') + ':' + String(e % 60).padStart(2, '0'); renderSession(); }
function status(t) { A.status = t; const s = $('#sess-status'); if (s) s.textContent = t; }
function show(id) { for (const s of ['home', 'session', 'results']) $('#' + s).classList.toggle('hidden', s !== id); scrollTo(0, 0); }

function renderSession() {
  const done = A.steps.filter((s) => s.st === 'done').length, fail = A.steps.some((s) => s.st === 'fail');
  $('#sess-step').textContent = t('step_n', Math.min(done + 1, 5));
  $('#sess-prog').style.width = (done / 5 * 100) + '%'; $('.progress').classList.toggle('fail', fail);
  $('#sess-title').textContent = t(fail ? 'st_attention' : done === 5 ? 'st_ready' : 'st_reducing');
  const g = { pending: '○', run: '◐', done: '●', fail: '✕' };
  $('#log').innerHTML = A.steps.map((s) => `<li class="${s.st}"><span class="t">${s.at || '--:--'}</span><span class="g">${g[s.st]}</span><span><b>${esc(t(s.k))}</b>${s.d ? `<small>${esc(s.d)}</small>` : ''}</span></li>`).join('');
  const act = $('#sess-actions'); act.innerHTML = '';
  if (fail) { const r = el('button', 'primary-btn', esc(t('retry'))); r.onclick = () => run(pipeline); const b = el('button', 'ghost', esc(t('back_home'))); b.onclick = home; act.append(r, b); }
  else { const c = el('button', 'ghost', esc(t('cancel'))); c.onclick = home; act.append(c); }
}
function home() { A.task++; clearSession(); show('home'); }
async function run(fn) { const id = ++A.task; show('session'); freshSteps(); try { await fn(() => id !== A.task); } catch (e) { if (id === A.task) { const i = A.steps.findIndex((s) => s.st === 'run'); step(i < 0 ? 0 : i, 'fail', String(e.message || e)); } } }

// ---------------------------------------------------------------- automatic pipeline
function fitsHeader(bytes) {
  const t = new TextDecoder('latin1').decode(bytes.slice(0, 2880 * 12)); if (!t.startsWith('SIMPLE')) return null; const h = {};
  for (let i = 0; i + 80 <= t.length; i += 80) { const c = t.slice(i, i + 80), k = c.slice(0, 8).trim(); if (k === 'END') break; if (c.slice(8, 10) !== '= ') continue;
    let v = c.slice(10); const q = v.indexOf("'"); v = q >= 0 ? v.slice(q + 1).split("'")[0] : v.split('/')[0]; h[k] = v.trim(); }
  return h;
}
function applyHeader() {
  const h = A.images[0] && fitsHeader(A.images[0].bytes); if (!h) return;
  if (!A.ra || !A.dec) {
    if (h.OBJCTRA && h.OBJCTDEC) { A.ra = h.OBJCTRA; A.dec = h.OBJCTDEC; }
    else if (isFinite(+h.RA) && isFinite(+h.DEC) && h.RA) { A.ra = fmt(+h.RA / 15); A.dec = fmt(+h.DEC, true); }
    else if (h.CRVAL1 && h.CRVAL2) { A.ra = fmt(+h.CRVAL1 / 15); A.dec = fmt(+h.CRVAL2, true); }
  }
  if (+h.FOCALLEN > 0) set('focal_length_mm', +h.FOCALLEN);
  if (+h.XPIXSZ > 0) { set('pixel_width_mm', +h.XPIXSZ / 1000); set('pixel_height_mm', (+(h.YPIXSZ || h.XPIXSZ)) / 1000); }
}
async function downloadCatalog(cancelled) {
  if (!A.ra || !A.dec) throw new Error(t('need_center'));
  const base = api.url(A.s, A.ra, A.dec);
  for (const h of VIZIER) { if (cancelled()) return;
    status(t('dl_gaia', h));
    try { const tx = await fetchRetry(base.replace('vizier.cds.unistra.fr', h), { tries: 2, min: 200, timeout: 40000 });
      if (tx.split('\n').filter((l) => /^[ 0-9]/.test(l)).length < 5) continue;
      A.M.FS.writeFile('/in/gaia.tsv', tx); A.catalog = '/in/gaia.tsv'; status(t('dl_gaia_ok', h, Math.round(tx.length / 1024))); return;
    } catch (e) { } }
  throw new Error(t('no_vizier'));
}
async function pipeline(cancelled) {
  if (!A.images.length) throw new Error(t('err_no_images'));
  step(0, 'done', t(A.images.length >= 3 ? 'n_images' : 'n_images_few', A.images.length));
  if (!A.ra || !A.dec) applyHeader();
  if (!A.ra || !A.dec) { step(1, 'fail', t('need_center')); openCenter(); return; }
  if (!A.catalog) { step(1, 'run'); await downloadCatalog(cancelled); }
  if (cancelled()) return; step(1, 'done', A.catalog.split('/').pop());
  step(2, 'run'); await new Promise((r) => setTimeout(r, 60)); if (cancelled()) return;
  const t0 = performance.now(), n = api.reduce(A.s, A.catalog, A.ra, A.dec, A.object);
  A.result = JSON.parse(api.json(A.s));
  if (n < 0 || A.result.error) throw new Error(coreError(A.result));
  const fs = frames(), rms = fs.reduce((a, f) => a + ((f.dra ?? 0) + (f.dde ?? 0)) / 2, 0) / fs.length;
  step(2, 'done', t('plate_done', Math.round(fs.reduce((a, f) => a + f.nref, 0) / fs.length), rms.toFixed(2), ((performance.now() - t0) / 1000).toFixed(2)));
  step(3, 'done', movers().length ? t('n_objects', movers().length) : t('none'));
  if (A.example) { step(4, 'done', t('example_offline')); showResults(); return; }
  step(4, 'run'); status(t('skybot_running')); await identify(); if (cancelled()) return;
  const lab = movers().map((m) => A.identified[m.desig]).filter(Boolean)[0];
  step(4, 'done', movers().length ? (lab ? lab.split('·')[0] : t('unidentified')) : t('n_known', A.known.length));
  showResults();
}
async function identify() {
  const f0 = frames()[0]; if (!f0) return; const code = get('mpc_code');
  const c = dbl2(api.center, 0);
  if (c) { try { A.known = await skybot(f0.jd, c[0], c[1], Math.max(f0.w, f0.h) * 1.1 / 3600 * .75, code); } catch (e) { A.known = []; } }
  for (const [i, m] of movers().entries()) {
    const q = m.meas[0];
    try { const l = await skybot(q.jd, q.rad, q.decd, 2 / 60, code); const a = l[0];
      if (!a) { A.identified[m.desig] = t('unidentified'); continue; }
      const d = A.object || a.number || a.name; api.desig(A.s, i, d);
      A.identified[d] = `${a.number ? `(${a.number}) ` : ''}${a.name} · ${a.cls} · V ${a.v.toFixed(1)} · ${a.sep.toFixed(0)}″`; m.desig = d;
    } catch (e) { A.identified[m.desig] = t('unidentified_off'); }
  }
}
async function startFiles(list) {
  clearSession(); for (const f of list) addFile(f.name, new Uint8Array(await f.arrayBuffer()));
  applyHeader(); run(pipeline);
}

// ---------------------------------------------------------------- results
const V = { k: 0, blink: true, layers: { track: true, known: true, refs: false, dets: false }, img: [], scale: 1, timer: 0 };
function showResults() {
  show('results'); V.k = 0; V.blink = true; V.img = frames().map((f, k) => frameImage(k));
  $('#frames').innerHTML = ''; frames().forEach((f, i) => { const b = el('button', '', String(i + 1)); b.onclick = () => { V.blink = false; V.k = i; draw(); }; $('#frames').append(b); });
  clearInterval(V.timer); V.timer = setInterval(() => { if (V.blink && !$('#results').classList.contains('hidden')) { V.k = (V.k + 1) % frames().length; draw(); } }, 600);
  renderSheet(); draw();
}
function frameImage(k) {
  const pw = A.M._malloc(8); api.size(A.s, k, pw, pw + 4); const w = A.M.getValue(pw, 'i32'), h = A.M.getValue(pw + 4, 'i32'); A.M._free(pw);
  const p = A.M._malloc(w * h * 4); api.render(A.s, k, p); const data = new Uint8ClampedArray(A.M.HEAPU8.buffer, p, w * h * 4).slice(); A.M._free(p);
  const c = document.createElement('canvas'); c.width = w; c.height = h; c.getContext('2d').putImageData(new ImageData(data, w, h), 0, 0); return c;
}
function draw() {
  const f = frames()[V.k]; if (!f) return; const cv = $('#sky'), wrap = $('#imgwrap');
  const base = Math.min(wrap.clientWidth / f.w, (wrap.clientHeight || innerHeight * .55) / f.h), s = base * V.scale, dpr = devicePixelRatio || 1;
  cv.style.width = f.w * s + 'px'; cv.style.height = f.h * s + 'px'; cv.width = f.w * s * dpr; cv.height = f.h * s * dpr;
  const g = cv.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0); g.imageSmoothingEnabled = false; g.drawImage(V.img[V.k], 0, 0, f.w * s, f.h * s);
  const P = (x, y) => [(x + .5) * s, (y + .5) * s];
  [...$('#frames').children].forEach((b, i) => b.classList.toggle('on', i === V.k));
  $('#ro-frame').textContent = t('plate_k', V.k + 1, f.time.slice(11));
  if (V.layers.dets) { g.strokeStyle = 'rgba(128,212,255,.35)'; g.lineWidth = .8; for (const d of f.dets) { const [x, y] = P(d[0], d[1]); g.beginPath(); g.arc(x, y, 3 * s, 0, 7); g.stroke(); } }
  if (V.layers.refs) { g.lineWidth = 1; for (const r of f.refs) { const [x, y] = P(r[0], r[1]); g.strokeStyle = r[2] ? 'rgba(155,229,100,.7)' : 'rgba(255,92,79,.7)'; reticle(g, x, y, 5 * s, .5); } }
  const idNames = movers().map((m) => A.identified[m.desig] || '');
  if (V.layers.known) { g.font = `500 ${Math.max(9, 3 * s)}px ui-monospace, Menlo, monospace`; g.textAlign = 'center';
    for (const a of A.known) { if (idNames.some((n) => n.includes(a.name))) continue; const q = dbl2(api.s2p, V.k, a.ra, a.dec); if (!q) continue;
      const [x, y] = P(q[0], q[1]); g.strokeStyle = 'rgba(128,212,255,.7)'; g.lineWidth = 1; reticle(g, x, y, 7 * s, .35); g.fillStyle = 'rgba(128,212,255,.85)'; g.fillText(a.name, x, y - 11 * s); } }
  for (const m of movers()) {
    const pts = m.meas.map((q) => P(q.x, q.y));
    if (V.layers.track) { g.setLineDash([3, 3]); g.strokeStyle = 'rgba(255,181,71,.5)'; g.lineWidth = 1.2; g.beginPath(); pts.forEach(([x, y], i) => i ? g.lineTo(x, y) : g.moveTo(x, y)); g.stroke(); g.setLineDash([]);
      g.fillStyle = 'rgba(255,181,71,.6)'; for (const [x, y] of pts) { g.beginPath(); g.arc(x, y, 1.5, 0, 7); g.fill(); } }
    if (pts[V.k]) { const [x, y] = pts[V.k]; g.strokeStyle = C.amber; g.lineWidth = 1.8; reticle(g, x, y, 10 * s);
      g.fillStyle = C.amber; g.font = `700 ${Math.max(11, 3.5 * s)}px ui-monospace, Menlo, monospace`; g.textAlign = 'center';
      g.fillText((A.identified[m.desig] || moverLabel(m)).split('·')[0].trim(), x, y + 20 * s); }
  }
  if (A.inspect) { const [x, y] = P(A.inspect.x, A.inspect.y); g.strokeStyle = '#fff'; g.lineWidth = 1.4; reticle(g, x, y, 8 * s, .6); }
}
function rows(r) { return `<div class="rows">${r.map(([k, v, c]) => `<div><span class="tag">${esc(k)}</span><span style="${c ? 'color:' + c : ''}">${esc(v)}</span></div>`).join('')}</div>`; }
function renderSheet() {
  const sh = $('#sheet'); let h = '';
  if (A.inspect) h += inspectHTML();
  if (A.example) h += `<section class="panel red"><span class="tag" style="color:var(--red)">${esc(t('example_banner'))}</span></section>`;
  if (!movers().length) h += `<section><span class="tag amber">${esc(t('result'))}</span><h2>${esc(t('nothing_moved'))}</h2><p class="status">${esc(t('nothing_moved_d', A.known.length))}</p></section>`;
  for (const m of movers()) {
    const id = A.identified[m.desig] || '', known = id && id !== t('unidentified') && id !== t('unidentified_off'), parts = id.split('·').map((x) => x.trim());
    const mag = m.meas.reduce((a, q) => a + q.mag, 0) / m.meas.length;
    h += `<section class="panel ${known ? 'amber' : 'red'}"><div class="bar"><span class="tag ${known ? 'amber' : ''}" style="${known ? '' : 'color:var(--red)'}">${esc(t('card_obj', known ? t('identified') : t('no_ident')))}</span><span style="color:${known ? C.amber : C.red}">${known ? '✓' : '?'}</span></div>
      <h2>${esc(known ? parts[0] : t('maybe_new'))}</h2>
      ${rows([...(known ? [[t('r_class'), parts[1]], [t('r_pred'), parts[2]], [t('r_dev'), parts[3]]] : []), ...(m.desig ? [] : [[t('r_temp'), moverLabel(m)]]), [t('r_rate'), (m.rate ?? 0).toFixed(2) + '″/min'], [t('r_pa'), (m.pa ?? 0).toFixed(1) + '°'], [t('r_mag'), mag.toFixed(1)]])}
      <p class="tag" style="margin-top:14px">${esc(t('measures'))}</p>
      <div class="meas">${m.meas.map((q, i) => `<div><b>${i + 1}</b><span class="d">${q.time.slice(11)}</span><span class="c">${q.ra}  ${q.dec}</span><span class="d">${q.mag.toFixed(1)}</span></div>`).join('')}</div>
      <div class="actions" style="margin-top:16px"><button class="primary-btn" data-dl="0">${esc(t('dl_ades'))}</button><button class="ghost" data-dl="all">${esc(t('dl_all'))}</button></div></section>`;
  }
  h += `<section class="panel"><span class="tag">${esc(t('plate_quality'))}</span>${frames().map((f, i) => { const r = ((f.dra ?? 0) + (f.dde ?? 0)) / 2;
    const [q, c] = f.nref < 6 ? [t('q_few'), C.red] : r < .25 ? [t('q_exc'), C.green] : r < .6 ? [t('q_good'), C.amber] : [t('q_weak'), C.red];
    return `<div class="q"><div class="h"><b>${esc(t('plate_n', i + 1))}</b><span style="color:${c}">${q}</span><span class="m">${f.nref} refs · ${r.toFixed(2)}″ · ${(f.seeing ?? 0).toFixed(1)}″</span></div><div class="bar2"><i style="width:${Math.max(5, Math.min(100, (1 - r) * 100))}%;background:${c}"></i></div></div>`; }).join('')}</section>`;
  if (A.known.length) h += `<section class="panel"><span class="tag">${esc(t('known_field'))}</span>${rows(A.known.slice(0, 8).map((a) => [(a.number ? `(${a.number}) ` : '') + a.name, 'V ' + a.v.toFixed(1), a.v < 18.5 ? C.cyan : 'var(--dim)']))}${A.known.length > 8 ? `<p class="tag">${esc(t('n_fainter', A.known.length - 8))}</p>` : ''}</section>`;
  h += `<button class="ghost" id="new-session">${esc(t('new_session'))}</button>`;
  sh.innerHTML = h;
  sh.querySelectorAll('[data-dl]').forEach((b) => b.onclick = () => download(b.dataset.dl));
  $('#new-session').onclick = home; const x = $('#insp-x'); if (x) x.onclick = () => { A.inspect = null; renderSheet(); draw(); };
}
function download(which) {
  const files = [['ADESReport.psv', 0], ['PhotReport.txt', 1], ['MPCReport.txt', 2]].filter((f, i) => which === 'all' || String(i) === which);
  for (const [name, k] of files) { const a = el('a'); a.href = URL.createObjectURL(new Blob([api.report(A.s, k)], { type: 'text/plain' })); a.download = name; a.click(); }
}
// "What is this?"
const simbadType = (k) => (SIMBAD_TYPES[k] ? SIMBAD_TYPES[k][LANG === 'pt' ? 1 : 0] : k);
function inspectHTML() {
  const i = A.inspect, coord = `${fmt(i.ra / 15)}   ${fmt(i.dec, true)}`;
  let body = `<p class="status">${esc(t('querying'))}</p>`;
  if (i.r) { const r = i.r; let title, det, col;
    if (r.ast) { title = (r.ast.number ? `(${r.ast.number}) ` : '') + r.ast.name; det = `${t('asteroid')} · ${r.ast.cls} · V ${r.ast.v.toFixed(1)}`; col = C.amber; }
    else if (r.vsx && r.vsx[0]) { const v = r.vsx[0]; title = v[0]; det = `${t('variable')} ${v[1]} · ${v[2]}–${v[3]} mag${v[4] ? ` · P = ${(+v[4]).toFixed(+v[4] < 10 ? 4 : 1)} d` : ''}`; col = C.red; }
    else if (r.sb && r.sb[0]) { title = r.sb[0].id; det = 'SIMBAD · ' + simbadType(r.sb[0].type); col = C.cyan; }
    else if (r.gaia) { title = t('gaia_star'); det = `G ${(+r.gaia[1]).toFixed(2)}${r.gaia[2] ? ` · BP−RP ${(+r.gaia[2]).toFixed(2)}` : ''}${+r.gaia[3] > 0 ? ' · ' + t('light_years', Math.round(3261.6 / +r.gaia[3])) : ''}`; col = C.green; }
    else { title = t('nothing_here'); det = t('nothing_here_d'); col = 'var(--dim)'; }
    const rr = []; if (r.ast) rr.push(['SkyBoT', `${r.ast.name} · ${r.ast.sep.toFixed(0)}″`]);
    (r.vsx || []).slice(0, 2).forEach((v) => rr.push(['VSX', `${v[0]} · ${v[1]}`])); (r.sb || []).slice(0, 2).forEach((x) => rr.push(['SIMBAD', `${x.id} · ${x.type} · ${x.sep.toFixed(1)}″`]));
    if (r.gaia) rr.push(['Gaia DR3', `G ${(+r.gaia[1]).toFixed(2)} · ${(+r.gaia[0]).toFixed(1)}″`]);
    body = `<div class="verdict">${esc(title)}</div><p class="status" style="color:${col}">${esc(det)}</p>${rows(rr)}${r.failed.length ? `<p class="tag">${esc(t('no_answer', r.failed.join(', ')))}</p>` : ''}`; }
  return `<section class="panel white"><div class="bar"><span class="tag" style="color:var(--text)">${esc(t('what_is'))}</span><button id="insp-x">✕</button></div><p class="tag">${coord}</p>${body}</section>`;
}
async function inspectAt(x, y) {
  const f = frames()[V.k]; let best = null, bd = 8;
  for (const d of f.dets) { const dd = Math.hypot(d[0] - x, d[1] - y); if (dd < bd) { bd = dd; best = d; } } if (best) { x = best[0]; y = best[1]; }
  const sky = dbl2(api.p2s, V.k, x, y); if (!sky) return;
  const id = Date.now(); A.inspect = { x, y, ra: sky[0], dec: sky[1], id }; V.blink = false; renderSheet(); draw();
  $('#ro-coord').textContent = `${fmt(sky[0] / 15)}  ${fmt(sky[1], true)}`;
  if (innerWidth <= 860) $('#sheet').scrollIntoView({ behavior: 'smooth' });
  const [sb, vx, ga, as] = await Promise.all([simbad(sky[0], sky[1], 6), vizier('B/vsx/vsx', sky[0], sky[1], 8, 'Name,Type,max,min,Period'),
    vizier('I/355/gaiadr3', sky[0], sky[1], 6, '_r,Gmag,BP-RP,Plx'), skybot(f.jd, sky[0], sky[1], 10 / 3600, get('mpc_code')).catch(() => null)]);
  if (!A.inspect || A.inspect.id !== id) return;
  const failed = []; if (!sb) failed.push('SIMBAD'); if (!vx) failed.push('VSX'); if (!ga) failed.push('Gaia'); if (!as) failed.push('SkyBoT');
  A.inspect.r = { sb: sb || [], vsx: (vx || []).filter((r) => r[0]), gaia: (ga || []).find((r) => isFinite(+r[1]) && r[1] !== ''), ast: as && as[0], failed };
  renderSheet(); draw();
}

// ---------------------------------------------------------------- dialogs
function openCenter() { const d = $('#center-dlg'); $('#c-ra').value = A.ra; $('#c-dec').value = A.dec; $('#c-msg').textContent = ''; d.showModal(); }
async function resolveName(n) {
  for (const b of ['https://vizier.cfa.harvard.edu/viz-bin/nph-sesame/-oI/A?', 'https://cds.unistra.fr/cgi-bin/nph-sesame/-oI/A?']) {
    try { const t = await fetchRetry(b + encodeURIComponent(n), { tries: 1, timeout: 20000 }); const l = t.split('\n').find((x) => x.startsWith('%J '));
      if (l) { const v = l.slice(3).trim().split(/\s+/).map(Number); return [v[0], v[1]]; } } catch (e) { } }
  return null;
}
function renderSettings() {
  const cur = A.profile;
  $('#profiles').innerHTML = profiles().map((p, i) => `<button type="button" class="profile" data-p="${i}"><span><b>${esc(p.name)}</b><small>${p.focal} mm · ${(p.pw * 1000).toFixed(1)}×${(p.ph * 1000).toFixed(1)} µm · ${(p.pw / p.focal * 206264.8).toFixed(2)}″/px · ${esc(p.mpc)}</small></span><span class="ck">${p.name === cur ? '✓' : ''}</span></button>`).join('');
  $('#profiles').querySelectorAll('[data-p]').forEach((b) => b.onclick = () => { applyProfile(profiles()[+b.dataset.p]); renderSettings(); });
  renderParams();
}

// parameter screen generated from astro_settings_schema(): collapsible groups, min/max validation
function renderParams() {
  const S = A.schema || { groups: [], settings: [] }, L = LANG === 'pt' ? 'pt' : 'en', box = $('#params'); box.innerHTML = '';
  for (const g of S.groups) {
    const fs = S.settings.filter((f) => f.group === g.id); if (!fs.length) continue;
    const d = el('details', 'group'); d.open = !g.collapsed; d.append(el('summary', '', esc(g[L]))); box.append(d);
    for (const f of fs) {
      const hint = f['hint_' + L], lab = el('label', 'param'); lab.title = hint;
      const name = el('span', 'pname', esc(f[L])); let inp;
      if (f.type === 'bool') { inp = el('input'); inp.type = 'checkbox'; inp.checked = Number(get(f.key) || f.default) !== 0; inp.onchange = () => set(f.key, inp.checked ? 1 : 0); }
      else if (f.type === 'enum') { inp = el('select'); for (const o of f.options) { const op = el('option', '', esc(o[L])); op.value = o.value; inp.append(op); }
        const cur = get(f.key) || f.default; inp.value = f.options.some((o) => o.value === cur) ? cur : (f.options.find((o) => Number(o.value) === Number(cur)) || {}).value ?? f.default; inp.onchange = () => set(f.key, inp.value); }
      else { inp = el('input'); inp.value = get(f.key); inp.placeholder = f.default;
        if (f.type === 'number' || f.type === 'int') { inp.inputMode = 'decimal';
          inp.oninput = () => { const v = inp.value.trim().replace(',', '.'), n = Number(v);
            const ok = v === '' || (Number.isFinite(n) && n >= f.min && n <= f.max && (f.type !== 'int' || Number.isInteger(n)));
            inp.classList.toggle('invalid', !ok); inp.title = ok ? hint : t('out_of_range', f.min, f.max); if (ok) set(f.key, v === '' ? f.default : v); }; }
        else inp.oninput = () => set(f.key, inp.value); }
      inp.dataset.k = f.key;
      const rs = el('button', 'reset', '↺'); rs.type = 'button'; rs.title = t('reset_field', f.default === '' ? '—' : f.default);
      rs.onclick = () => { set(f.key, f.default); renderSettings(); };
      lab.append(name, inp, rs); lab.append(el('small', 'hint', esc(hint))); d.append(lab);
    }
  }
}

// ---------------------------------------------------------------- built-in example (offline)
function coreError(r) {
  const m = { no_images: 'err_no_images', image_unreadable: 'err_image_unreadable', catalog_empty: 'err_catalog_empty', internal: 'err_internal' };
  return r && m[r.code] ? t(m[r.code], r.detail || '') : (r && r.error) || t('reduce_error');
}
async function startExample() {
  clearSession(); show('session'); freshSteps(); status(t('example_loading'));
  try {
    const get1 = async (n, txt) => { const r = await fetch('example/' + n); if (!r.ok) throw new Error(n + ': HTTP ' + r.status); return txt ? r.text() : new Uint8Array(await r.arrayBuffer()); };
    const cfg = await get1('example.cfg', true); A.M.FS.writeFile('/in/example.cfg', cfg); api.loadCfg(A.s, '/in/example.cfg');
    const imgs = get('example_images').split(';').filter(Boolean), cat = get('example_catalog');
    for (const n of imgs) addFile(n, await get1(n));
    A.M.FS.writeFile('/in/' + cat, await get1(cat, true)); A.catalog = '/in/' + cat;
    A.ra = get('example_ra'); A.dec = get('example_dec');
  } catch (e) { step(0, 'fail', t('example_fail', e.message || e)); return; }
  A.example = get('example_no_submit') !== '0';
  run(pipeline);
}

// ---------------------------------------------------------------- events
addEventListener('DOMContentLoaded', () => {
  drawChart($('#chart')); addEventListener('resize', () => { drawChart($('#chart')); if (!$('#results').classList.contains('hidden')) draw(); });
  // animated sample plate
  const sc = $('#sample-anim'); let ph = 0;
  setInterval(() => { drawChart(sc, 11, false); const g = sc.getContext('2d'), w = sc.clientWidth, h = sc.clientHeight; ph ^= 1;
    const x = w * (ph ? .58 : .55), y = h * (ph ? .43 : .47); g.setLineDash([3, 4]); g.strokeStyle = 'rgba(255,181,71,.4)'; g.beginPath(); g.moveTo(w * .49, h * .55); g.lineTo(w * .64, h * .35); g.stroke(); g.setLineDash([]);
    g.strokeStyle = C.amber; g.lineWidth = 1.5; reticle(g, x, y, 15); g.fillStyle = '#fff'; g.beginPath(); g.arc(x, y, 1.6, 0, 7); g.fill(); }, 700);
  applyI18n();
  $('#lang').value = LANG; $('#lang').onchange = () => { setLang($('#lang').value); renderSettings(); if (A.steps.length) renderSession(); if (!$('#results').classList.contains('hidden')) { renderSheet(); draw(); } const b = $('#btn-settings'); if (b) b.textContent = profileName() + ' ⚙︎'; };
  $('#p-defaults').onclick = () => { for (const f of A.schema.settings) set(f.key, f.default); renderSettings(); };
  $('#go-files').onclick = () => $('#files').click();
  $('#go-example').onclick = () => startExample();
  $('#files').onchange = (e) => { if (e.target.files.length) startFiles([...e.target.files]); e.target.value = ''; };
  $('#sess-close').onclick = home; $('#res-back').onclick = () => { clearInterval(V.timer); home(); };
  $('#res-blink').onclick = () => { V.blink = !V.blink; $('#res-blink').textContent = V.blink ? '❚❚' : '▶'; };
  $('#res-layers').onclick = () => $('#layers').classList.toggle('hidden');
  $('#layers').querySelectorAll('input').forEach((i) => i.onchange = () => { V.layers[i.dataset.l] = i.checked; draw(); });
  $('#sky').onclick = (ev) => { const f = frames()[V.k], r = ev.target.getBoundingClientRect(), s = r.width / f.w; inspectAt((ev.clientX - r.left) / s - .5, (ev.clientY - r.top) / s - .5); };
  $('#imgwrap').addEventListener('wheel', (ev) => { if (!ev.ctrlKey && !ev.metaKey) return; ev.preventDefault(); V.scale = Math.max(1, Math.min(6, V.scale * (ev.deltaY < 0 ? 1.15 : .87))); draw(); }, { passive: false });
  $('#sky').ondblclick = () => { V.scale = V.scale > 1 ? 1 : 2.5; draw(); };
  // drag and drop files onto the page
  addEventListener('dragover', (e) => e.preventDefault());
  addEventListener('drop', (e) => { e.preventDefault(); if (e.dataTransfer.files.length) startFiles([...e.dataTransfer.files]); });
  // dialogs
  $('#btn-settings').onclick = () => { renderSettings(); $('#settings-dlg').showModal(); };
  $('#p-save').onclick = () => { const n = $('#p-name').value.trim(); if (!n) return;
    const p = { name: n, focal: +get('focal_length_mm'), pw: +get('pixel_width_mm'), ph: +get('pixel_height_mm'), flip: +get('flip_horizontal'), mpc: get('mpc_code') };
    const user = profiles().filter((x) => x.name !== n); user.push(p); lsSet('profiles', JSON.stringify(user));
    applyProfile(p); $('#p-name').value = ''; renderSettings(); };
  $('#c-find').onclick = async () => { $('#c-msg').textContent = t('searching'); const r = await resolveName($('#c-name').value);
    if (r) { $('#c-ra').value = fmt(r[0] / 15); $('#c-dec').value = fmt(r[1], true); $('#c-msg').textContent = t('found'); } else $('#c-msg').textContent = t('not_found'); };
  $('#c-cancel').onclick = () => $('#center-dlg').close();
  $('#c-ok').onclick = () => { A.ra = $('#c-ra').value.trim(); A.dec = $('#c-dec').value.trim(); if (!A.ra || !A.dec) return; $('#center-dlg').close(); run(pipeline); };
  boot().catch((e) => { document.body.insertAdjacentHTML('beforeend', `<p class="status" style="position:fixed;bottom:10px;left:10px;z-index:9;color:var(--red)">${esc(t('core_fail', e.message || e))}</p>`); });
});

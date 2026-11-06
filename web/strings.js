// Piazzi web — todas as strings da interface (inglês padrão, português opcional).
// t('chave', a, b…) troca {0}, {1}… pelos argumentos.
'use strict';
const STRINGS = {
  // início
  tagline_tag: ['Pocket observatory', 'Observatório de bolso'],
  lede: ['Measure the sky. Find what moves.', 'Meça o céu. Encontre o que se move.'],
  card_files_tag: ['Your plates', 'Suas placas'],
  card_files_title: ['Open my images', 'Abrir minhas imagens'],
  card_files_desc: ['FITS or SBIG from your telescope. Catalog, plate, motion and identification — automatic.', 'FITS ou SBIG do seu telescópio. Catálogo, placa, movimento e identificação — automático.'],
  card_example_tag: ['Try it', 'Experimente'],
  card_example_title: ['Open example', 'Abrir exemplo'],
  card_example_desc: ['Three ZTF images of asteroid (5157) Hindemith with a Gaia DR3 catalog — runs fully offline.', 'Três imagens do ZTF do asteroide (5157) Hindemith com catálogo Gaia DR3 — roda totalmente offline.'],
  default_profile: ['Default', 'Padrão'],
  core_fail: ['Failed to load the core: {0}', 'Falha ao carregar o núcleo: {0}'],
  // sessão
  close: ['Close', 'Fechar'],
  step_n: ['Step {0}/5', 'Etapa {0}/5'],
  log_title: ['Observing log', 'Diário de observação'],
  st_attention: ['Needs attention.', 'Precisa de atenção.'],
  st_ready: ['Done.', 'Pronto.'],
  st_reducing: ['Reducing the plates…', 'Reduzindo as placas…'],
  retry: ['Try again', 'Tentar de novo'],
  back_home: ['Back to start', 'Voltar ao início'],
  cancel: ['Cancel', 'Cancelar'],
  s_images: ['Images', 'Imagens'],
  s_catalog: ['Gaia DR3 catalog', 'Catálogo Gaia DR3'],
  s_plate: ['Astrometric plate', 'Placa astrométrica'],
  s_motion: ['Moving objects', 'Objetos em movimento'],
  s_ident: ['Identification', 'Identificação'],
  n_images: ['{0} images', '{0} imagens'],
  n_images_few: ['{0} image(s) — motion needs 3+', '{0} imagem(ns) — movimento precisa de 3+'],
  need_center: ['Enter the field centre', 'Informe o centro do campo'],
  dl_gaia: ['Downloading Gaia DR3 from {0}…', 'Baixando Gaia DR3 de {0}…'],
  dl_gaia_ok: ['Gaia DR3 catalog downloaded from {0} ({1} kB).', 'Catálogo Gaia DR3 baixado de {0} ({1} kB).'],
  no_vizier: ['No VizieR mirror answered — check the internet connection or add a Gaia DR3 TSV with the images.', 'Nenhum espelho do VizieR respondeu — verifique a internet ou junte um TSV do Gaia DR3 às imagens.'],
  plate_done: ['{0} refs per image · {1}″ · {2} s', '{0} refs por imagem · {1}″ · {2} s'],
  n_objects: ['{0} object(s)', '{0} objeto(s)'],
  none: ['none', 'nenhum'],
  skybot_running: ['Identifying with SkyBoT…', 'Identificando no SkyBoT…'],
  unidentified: ['not identified', 'não identificado'],
  unidentified_off: ['not identified (SkyBoT unavailable)', 'não identificado (SkyBoT indisponível)'],
  n_known: ['{0} known asteroids in the field', '{0} asteroides conhecidos no campo'],
  example_offline: ['example: offline, SkyBoT skipped', 'exemplo: offline, SkyBoT não consultado'],
  example_loading: ['Loading the example…', 'Carregando o exemplo…'],
  example_fail: ['Could not load the example files: {0}', 'Não consegui carregar os arquivos do exemplo: {0}'],
  reduce_error: ['Reduction error', 'Erro na redução'],
  // erros do núcleo (código → mensagem)
  err_no_images: ['No images — add at least two exposures of the same field.', 'Nenhuma imagem — adicione pelo menos duas exposições do mesmo campo.'],
  err_image_unreadable: ['Image unreadable: {0} — use FITS or SBIG files.', 'Imagem ilegível: {0} — use arquivos FITS ou SBIG.'],
  err_catalog_empty: ['Catalog empty or unreadable — check the internet connection or load a Gaia DR3 TSV.', 'Catálogo vazio ou ilegível — verifique a internet ou carregue um TSV do Gaia DR3.'],
  err_internal: ['Internal error: {0}', 'Erro interno: {0}'],
  // resultado
  back: ['Back', 'Voltar'],
  blink: ['Blink', 'Piscar'],
  layers: ['Layers', 'Camadas'],
  l_track: ['Object track', 'Rastro do objeto'],
  l_known: ['Known asteroids', 'Asteroides conhecidos'],
  l_refs: ['Reference stars', 'Estrelas de referência'],
  l_dets: ['All detections', 'Todas as detecções'],
  tap_identify: ['tap = identify', 'toque = identificar'],
  plate_k: ['Plate {0} · {1}', 'Placa {0} · {1}'],
  object_n: ['object {0}', 'objeto {0}'],
  result: ['Result', 'Resultado'],
  nothing_moved: ['Nothing moved.', 'Nada se moveu.'],
  nothing_moved_d: ['The plates were solved, but no object moved between the images above the limit. {0} known asteroid(s) cross the field, probably too faint.', 'As placas foram resolvidas, mas nenhum objeto andou entre as imagens acima do limite. {0} asteroide(s) conhecido(s) passam no campo, provavelmente fracos demais.'],
  card_obj: ['Object card · {0}', 'Ficha do objeto · {0}'],
  identified: ['identified', 'identificado'],
  no_ident: ['no identification', 'sem identificação'],
  maybe_new: ['Possible new object', 'Possível objeto novo'],
  r_class: ['Class', 'Classe'],
  r_pred: ['Predicted brightness', 'Brilho previsto'],
  r_dev: ['Offset from ephemeris', 'Desvio da efeméride'],
  r_temp: ['Temporary designation', 'Designação temporária'],
  r_rate: ['Rate', 'Velocidade'],
  r_pa: ['Position angle', 'Ângulo de posição'],
  r_mag: ['Measured magnitude', 'Magnitude medida'],
  measures: ['Measurements', 'Medidas'],
  dl_ades: ['Download ADES report (MPC)', 'Baixar relatório ADES (MPC)'],
  dl_all: ['ADES · photometric · MPC 80', 'ADES · fotométrico · MPC 80'],
  example_banner: ['Example data — not your observation. Do not submit these reports to the MPC.', 'Dados de exemplo — não é observação sua. Não envie estes relatórios ao MPC.'],
  plate_quality: ['Plate quality', 'Qualidade das placas'],
  q_few: ['few stars', 'poucas estrelas'], q_exc: ['excellent', 'excelente'], q_good: ['good', 'boa'], q_weak: ['weak', 'fraca'],
  plate_n: ['Plate {0}', 'Placa {0}'],
  known_field: ['Known asteroids in the field · SkyBoT', 'Asteroides conhecidos no campo · SkyBoT'],
  n_fainter: ['+ {0} fainter', '+ {0} mais fracos'],
  new_session: ['New session', 'Nova sessão'],
  // o que é isso?
  what_is: ['What is this?', 'O que é isso?'],
  querying: ['Querying SIMBAD, VSX, Gaia and SkyBoT…', 'Consultando SIMBAD, VSX, Gaia e SkyBoT…'],
  asteroid: ['Asteroid', 'Asteroide'],
  variable: ['Variable star', 'Estrela variável'],
  gaia_star: ['Gaia DR3 star', 'Estrela Gaia DR3'],
  light_years: ['~{0} light-years', '~{0} anos-luz'],
  nothing_here: ['Nothing known here', 'Nada conhecido aqui'],
  nothing_here_d: ['No catalogued object within 6″ — noise, artefact or something new.', 'Nenhum objeto catalogado num raio de 6″ — ruído, artefato ou algo novo.'],
  no_answer: ['no answer: {0}', 'sem resposta: {0}'],
  // centro do campo
  center_tag: ['Field centre', 'Centro do campo'],
  center_q: ['Where was the telescope pointing?', 'Para onde o telescópio apontava?'],
  center_d: ['The coordinates are not in the header. Search for the target name or type RA and Dec.', 'Não encontrei as coordenadas no cabeçalho. Procure pelo nome do alvo ou digite AR e Dec.'],
  ra_hms: ['RA h m s', 'AR h m s'],
  continue: ['Continue', 'Continuar'],
  searching: ['searching…', 'procurando…'],
  found: ['found', 'encontrado'],
  not_found: ['not found', 'não encontrado'],
  // ajustes
  settings: ['Settings', 'Ajustes'],
  language: ['Language', 'Idioma'],
  equipment: ['Equipment', 'Equipamento'],
  profile_ph: ['Profile name (e.g. C11 + ASI533)', 'Nome do perfil (ex.: C11 + ASI533)'],
  parameters: ['Parameters', 'Parâmetros'],
  restore_defaults: ['Restore defaults', 'Restaurar padrões'],
  reset_field: ['Restore default ({0})', 'Restaurar padrão ({0})'],
  out_of_range: ['Must be between {0} and {1}', 'Deve estar entre {0} e {1}'],
};
const SIMBAD_TYPES = {
  '*': ['Star', 'Estrela'], 's*r': ['Red supergiant', 'Supergigante vermelha'], G: ['Galaxy', 'Galáxia'], QSO: ['Quasar', 'Quasar'],
  PN: ['Planetary nebula', 'Nebulosa planetária'], 'V*': ['Variable star', 'Estrela variável'], 'EB*': ['Eclipsing binary', 'Binária eclipsante'],
  '**': ['Double star', 'Estrela dupla'], 'PM*': ['High proper-motion star', 'Estrela de alto movimento próprio'], HII: ['HII region', 'Região HII'],
  'Cl*': ['Cluster', 'Aglomerado'],
};
let LANG = 'en';
try { LANG = localStorage.getItem('lang') || ((navigator.language || '').toLowerCase().startsWith('pt') ? 'pt' : 'en'); } catch (e) { LANG = 'en'; }
if (LANG !== 'pt') LANG = 'en';
function t(k, ...a) {
  const e = STRINGS[k]; let s = e ? e[LANG === 'pt' ? 1 : 0] : k;
  a.forEach((v, i) => { s = s.split('{' + i + '}').join(String(v)); });
  return s;
}
function setLang(l) { LANG = l === 'pt' ? 'pt' : 'en'; try { localStorage.setItem('lang', LANG); } catch (e) { /* sem armazenamento */ } applyI18n(); }
// textos estáticos do HTML: data-i18n (texto), data-i18n-ph (placeholder), data-i18n-aria (aria-label)
function applyI18n() {
  document.documentElement.lang = LANG === 'pt' ? 'pt-BR' : 'en';
  document.querySelectorAll('[data-i18n]').forEach((e) => { e.textContent = t(e.dataset.i18n); });
  document.querySelectorAll('[data-i18n-ph]').forEach((e) => { e.placeholder = t(e.dataset.i18nPh); });
  document.querySelectorAll('[data-i18n-aria]').forEach((e) => { e.setAttribute('aria-label', t(e.dataset.i18nAria)); });
}

// Tabela de textos en -> pt do app de desktop e escolha de idioma (guardada nas preferências).
#include "strings.hpp"
#include "app.hpp"
#include <fstream>
#include <string>
#include <unordered_map>

namespace {
std::string g_lang = "en";
const std::unordered_map<std::string, const char *> &table() {
    static const std::unordered_map<std::string, const char *> T = {
        {"Images", "Imagens"},
        {"Reference catalog", "Catálogo de referência"},
        {"Astrometric plate", "Placa astrométrica"},
        {"Moving objects", "Objetos em movimento"},
        {"Identification", "Identificação"},
        {"local orbit", "órbita local"},
        {"%s · local orbit · %.0f\"", "%s · órbita local · %.0f\""},
        {"SkyBoT did not answer (no internet?) — using the local orbits.",
         "SkyBoT não respondeu (sem internet?) — usando as órbitas locais."},
        {"SkyBoT did not answer (no internet?) — no identification.", "SkyBoT não respondeu (sem internet?) — sem identificação."},
        {"No plate solved — SkyBoT not queried.", "Nenhuma placa resolvida — sem consulta ao SkyBoT."},
        {"No images — add at least two exposures of the same field", "Nenhuma imagem"},
        {" images", " imagens"},
        {" image(s) — motion needs 3+", " imagem(ns) — movimento precisa de 3+"},
        {"Enter the field centre", "Informe o centro do campo"},
        {"Reducing the plates…", "Reduzindo as placas…"},
        {"%d refs per image · %.2f\" · %.2f s", "%d refs por imagem · %.2f\" · %.2f s"},
        {"none", "nenhum"},
        {" object(s)", " objeto(s)"},
        {"Identifying…", "Identificando…"},
        {" known asteroids in the field", " asteroides conhecidos no campo"},
        {"Error: ", "Erro: "},
        {"Unexpected error in the worker thread.", "Erro inesperado na thread de trabalho."},
        {"Your plates", "Suas placas"},
        {"Open my images", "Abrir minhas imagens"},
        {"FITS or SBIG from your telescope. Or drag the files onto the window.",
         "FITS ou SBIG do seu telescópio. Ou arraste os arquivos para a janela."},
        {"Choose the images (and, optionally, the .cfg and the catalog)", "Escolha as imagens (e, se quiser, o .cfg e o catálogo)"},
        {"Pocket observatory", "Observatório de bolso"},
        {"  ·  settings", "  ·  ajustes"},
        {"Measure the sky. Find what moves.", "Meça o céu. Encontre o que se move."},
        {"Step ", "Etapa "},
        {"Needs attention.", "Precisa de atenção."},
        {"Done.", "Pronto."},
        {"Observing log", "Diário de observação"},
        {"TRY AGAIN", "TENTAR DE NOVO"},
        {"BACK TO START", "VOLTAR AO INÍCIO"},
        {"CANCEL", "CANCELAR"},
        {"↑ up", "↑ acima"},
        {"Ctrl/⌘ + click selects several · double-click opens a folder", "Ctrl/⌘ + clique seleciona vários · duplo clique abre pasta"},
        {"double-click opens a folder", "duplo clique abre pasta"},
        {"Cancel", "Cancelar"},
        {"Use folder", "Usar pasta"},
        {"Open", "Abrir"},
        {"mode", "modo"},
        {"Stack", "Empilhar"},
        {"Stacked: plate ", "Empilhado: placa "},
        {"Measurements on the stacked plate get note K.", "Medidas na placa empilhada recebem nota K."},
        {"black %", "preto %"},
        {"white %", "branco %"},
        {"stretch", "escala"},
        {"invert", "inverter"},
        {"R reference · S star · T trail fragment · L trail · ? loose",
         "R referência · S estrela · T fragmento de rastro · L rastro · ? solta"},
        {"header", "cabeçalho"},
        {"correct (s)", "corrigir (s)"},
        {"apply", "aplicar"},
        {"Time changed — reduce again to apply.", "Horário alterado — reduza de novo para aplicar."},
        {"Reduce again", "Reduzir de novo"},
        {"Save FITS with WCS (current plate)", "Salvar FITS com WCS (placa atual)"},
        {"image", "imagem"},
        {"FITS saved to ", "FITS salvo em "},
        {"Could not save the FITS.", "Falha ao salvar o FITS."},
        {"Contrast", "Contraste"},
        {"Detection list (Image Catalog)", "Lista de detecções (Image Catalog)"},
        {"Header and time", "Cabeçalho e horário"},
        {"Save image", "Salvar imagem"},
        {"pause", "pausar"},
        {"blink", "piscar"},
        {"layers", "camadas"},
        {"Object track", "Rastro do objeto"},
        {"Known asteroids", "Asteroides conhecidos"},
        {"Reference stars", "Estrelas de referência"},
        {"All detections", "Todas as detecções"},
        {"Plate ", "Placa "},
        {"new object?", "objeto novo?"},
        {"Reports saved to ", "Relatórios gravados em "},
        {"What is this?", "O que é isso?"},
        {"Querying SIMBAD, VSX, Gaia and SkyBoT…", "Consultando SIMBAD, VSX, Gaia e SkyBoT…"},
        {"Manual measurement · plate ", "Medida manual · placa "},
        {"RA / Dec", "AR / Dec"},
        {"SNR peak / total", "SNR pico / total"},
        {"CREATE OBJECT FROM THIS MEASUREMENT", "CRIAR OBJETO COM ESTA MEDIDA"},
        {"ADD TO THE MANUAL OBJECT", "ACRESCENTAR AO OBJETO MANUAL"},
        {"Measurement added to the manual object.", "Medida acrescentada ao objeto manual."},
        {"Object sheet · identified", "Ficha do objeto · identificado"},
        {"Object sheet · not identified", "Ficha do objeto · sem identificação"},
        {"Possible new object", "Possível objeto novo"},
        {"object %zu", "objeto %zu"},
        {"Label", "Rótulo"},
        {"Class", "Classe"},
        {"Detail", "Detalhe"},
        {"Rate", "Velocidade"},
        {"Position angle", "Ângulo de posição"},
        {"Measurements · MPC note · exclude", "Medidas · nota MPC · excluir"},
        {"SUBMIT TO MPC", "ENVIAR AO MPC"},
        {"SAVE ADES · PHOTOMETRIC · MPC 80", "SALVAR ADES · FOTOMÉTRICO · MPC 80"},
        {"Plate quality", "Qualidade das placas"},
        {"few stars", "poucas estrelas"},
        {"excellent", "excelente"},
        {"good", "boa"},
        {"poor", "fraca"},
        {"Plate %zu%s", "Placa %zu%s"},
        {" (stacked)", " (empilhada)"},
        {"Result", "Resultado"},
        {"Nothing moved.", "Nada se moveu."},
        {"No object moved between the images above the threshold. %d known asteroid(s) in the field.",
         "Nenhum objeto andou entre as imagens acima do limite. %d asteroide(s) conhecido(s) no campo."},
        {"Known asteroids in the field", "Asteroides conhecidos no campo"},
        {"NEW SESSION", "NOVA SESSÃO"},
        {"Downloading catalog from ", "Baixando catálogo de "},
        {"Catalog downloaded from ", "Catálogo baixado de "},
        {"No VizieR mirror answered (no internet?). Load a catalog file (Gaia DR3 TSV).",
         "Nenhum espelho do VizieR respondeu (sem internet?). Use um arquivo de catálogo."},
        {"Asteroid · ", "Asteroide · "},
        {"Variable star ", "Estrela variável "},
        {"Gaia DR3 star", "Estrela Gaia DR3"},
        {"No network answer", "Sem resposta da rede"},
        {"SkyBoT, VizieR and SIMBAD did not answer (no internet?).", "SkyBoT, VizieR e SIMBAD não responderam (sem internet?)."},
        {"Nothing known here", "Nada conhecido aqui"},
        {"No catalogued object within 6\".", "Nenhum objeto catalogado num raio de 6\"."},
        // ---------------------------------------------------------------- início, exemplo e idioma
        {"Example", "Exemplo"},
        {"Open example", "Abrir exemplo"},
        {"Asteroid (5157) Hindemith in three ZTF images. Runs offline.",
         "O asteroide (5157) Hindemith em três imagens do ZTF. Roda offline."},
        {"Example files not found: ", "Arquivos do exemplo não encontrados: "},
        {"Default", "Padrão"},
        {"Language", "Idioma"},
        // ---------------------------------------------------------------- centro do campo
        {"Field centre", "Centro do campo"},
        {"Where was the telescope pointing?", "Para onde o telescópio apontava?"},
        {"No coordinates in the header. Search for the target by name or type RA and Dec.",
         "Não encontrei as coordenadas no cabeçalho. Procure pelo nome do alvo ou digite AR e Dec."},
        {"found", "encontrado"},
        {"not found", "não encontrado"},
        {"no answer from Sesame (no internet?)", "sem resposta do Sesame (sem internet?)"},
        {"Searching…", "Buscando…"},
        {"Search", "Buscar"},
        {"RA h m s", "AR h m s"},
        {"CONTINUE", "CONTINUAR"},
        // ---------------------------------------------------------------- ajustes
        {"Settings", "Ajustes"},
        {"Equipment", "Equipamento"},
        {"Profile name (e.g. C11 + ASI533)", "Nome do perfil (ex.: C11 + ASI533)"},
        {"Save current", "Salvar atual"},
        {"Invalid equipment parameters — profile not saved.", "Parâmetros do equipamento inválidos — perfil não salvo."},
        {"Parameters", "Parâmetros"},
        {"Restore defaults", "Restaurar padrões"},
        {"Restore default", "Restaurar padrão"},
        {"default", "padrão"},
        {"out of range", "fora do intervalo"},
        {"not a number", "não é um número"},
        {"Defaults restored.", "Padrões restaurados."},
        {"(none)", "(nenhum)"},
        {"choose", "escolher"},
        {"remove", "remover"},
        {"Wait for the session to finish before changing this file.", "Espere a sessão terminar para trocar este arquivo."},
        {"Loading orbits…", "Carregando órbitas…"},
        {" orbits loaded.", " órbitas carregadas."},
        {"Session running — settings are locked until it finishes.", "Sessão em andamento — os ajustes ficam bloqueados até ela terminar."},
        {"Calibration and orbits", "Calibração e órbitas"},
        {"%d orbits loaded (offline identification).", "%d órbitas carregadas (identificação offline)."},
        {"Output", "Saída"},
        {"choose folder", "escolher pasta"},
        {"Folder for reports and FITS", "Pasta para relatórios e FITS"},
        // ---------------------------------------------------------------- listas (combos)
        {"sum", "soma"},
        {"mean", "média"},
        {"median", "mediana"},
        {"square root", "raiz"},
        {"click: identify", "clique: identificar"},
        {"click: measure", "clique: medir"},
        // ---------------------------------------------------------------- envio ao MPC
        {"Submission to the Minor Planet Center", "Envio ao Minor Planet Center"},
        {"Submit the ADES report?", "Enviar o relatório ADES?"},
        {"The file ADESReport.psv will be sent to the MPC through the submit_psv interface. Review the measurements first. In test "
         "mode the MPC only validates and answers, without recording.",
         "O arquivo ADESReport.psv será enviado ao MPC pela interface submit_psv. Revise as medidas antes. No modo teste o MPC "
         "só valida e responde, sem registrar."},
        {"e-mail for confirmation", "e-mail para confirmação"},
        {"you@example.com", "voce@exemplo.com"},
        {"invalid e-mail", "e-mail inválido"},
        {"Test mode (recommended the first time)", "Modo teste (recomendado na primeira vez)"},
        {"Close", "Fechar"},
        {"SUBMIT (TEST)", "ENVIAR (TESTE)"},
        {"SUBMIT FOR REAL", "ENVIAR DE VERDADE"},
        {"MPC answer: ", "Resposta do MPC: "},
        {"Submission failed (no internet?).", "Falha no envio (sem internet?)."},
        {"Example data — not your observation", "Dados de exemplo — não são uma observação sua"},
        {"These measurements come from the bundled ZTF example and must not be sent to the MPC. Open your own images to submit.",
         "Estas medidas vêm do exemplo do ZTF que acompanha o app e não podem ser enviadas ao MPC. Abra suas próprias imagens para "
         "enviar."},
        // ---------------------------------------------------------------- erros do núcleo
        {"No images — add at least two exposures of the same field", "Nenhuma imagem — acrescente ao menos duas exposições do mesmo campo"},
        {"Image unreadable: ", "Imagem ilegível: "},
        {" — use FITS or SBIG files", " — use arquivos FITS ou SBIG"},
        {"Catalog empty or unreadable — check the internet connection or load a Gaia DR3 TSV",
         "Catálogo vazio ou ilegível — verifique a conexão com a internet ou carregue um TSV do Gaia DR3"},
        {"Internal error: ", "Erro interno: "},
    };
    return T;
}
// valor de uma chave de string no JSON simples da API C
std::string jstr(const std::string &j, const std::string &k) {
    auto p = j.find("\"" + k + "\":\"");
    if (p == std::string::npos)
        return "";
    p += k.size() + 4;
    std::string out;
    for (size_t i = p; i < j.size() && j[i] != '"'; ++i) {
        if (j[i] == '\\' && i + 1 < j.size())
            ++i;
        out += j[i];
    }
    return out;
}
} // namespace

const char *tr(const char *en) {
    if (g_lang != "pt" || !en)
        return en;
    auto it = table().find(en);
    return it == table().end() ? en : it->second;
}
const std::string &lang() {
    return g_lang;
}
bool lang_pt() {
    return g_lang == "pt";
}
void set_lang(const std::string &l) {
    g_lang = l == "pt" ? "pt" : "en";
    std::ofstream(prefs_dir() + "lang.txt") << g_lang << "\n";
}
void load_lang() {
    std::ifstream f(prefs_dir() + "lang.txt");
    std::string l;
    if (f >> l)
        g_lang = l == "pt" ? "pt" : "en";
}

std::string core_error_message(const std::string &j) {
    std::string code = jstr(j, "code"), detail = jstr(j, "detail");
    if (code == "no_images")
        return tr("No images — add at least two exposures of the same field");
    if (code == "image_unreadable")
        return std::string(tr("Image unreadable: ")) + detail + tr(" — use FITS or SBIG files");
    if (code == "catalog_empty")
        return tr("Catalog empty or unreadable — check the internet connection or load a Gaia DR3 TSV");
    if (code == "internal")
        return std::string(tr("Internal error: ")) + detail;
    std::string e = jstr(j, "error");
    return e.empty() ? j.substr(0, 160) : e;
}

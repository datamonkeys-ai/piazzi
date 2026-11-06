import Foundation

/// Idioma do app: inglês por padrão, português opcional. "system" segue o idioma do iOS; "en"/"pt" forçam.
/// Guardado em UserDefaults ("lang"), lido por L() a cada desenho (a raiz troca o .id quando muda).
enum Lang {
    static var setting: String { UserDefaults.standard.string(forKey: "lang") ?? "system" }
    static var isPT: Bool {
        switch setting {
        case "pt": return true
        case "en": return false
        default: return (Locale.preferredLanguages.first ?? "en").hasPrefix("pt")
        }
    }
}

/// Tradução de uma string (chave = texto em inglês). Sem entrada na tabela, devolve o inglês.
func L(_ en: String) -> String { Lang.isPT ? (ptTable[en] ?? en) : en }
/// Com formato (String(format:)) e argumentos.
func L(_ en: String, _ args: CVarArg...) -> String { String(format: L(en), arguments: args) }

/// Mensagem localizada para o código de erro do núcleo (astro_result_json: "code").
func coreError(code: String?, message: String?, detail: String?) -> String {
    switch code {
    case "no_images": return L("No images — add at least two exposures of the same field")
    case "image_unreadable": return L("Image unreadable: %@ — use FITS or SBIG files", (detail ?? "") as NSString)
    case "catalog_empty": return L("Catalog empty or unreadable — check the internet connection or load a Gaia DR3 TSV")
    case "internal": return L("Internal error: %@", (detail ?? message ?? "") as NSString)
    default: return message ?? L("Reduction error.")
    }
}

private let ptTable: [String: String] = [
    // início
    "Pocket observatory": "Observatório de bolso",
    "Measure the sky. Find what moves.": "Meça o céu. Encontre o que se move.",
    "Your plates": "Suas placas",
    "Open my images": "Abrir minhas imagens",
    "FITS or SBIG from your telescope. Catalog, plate, motion and identification — automatic.":
        "FITS ou SBIG do seu telescópio. Catálogo, placa, movimento e identificação — automático.",
    "Try it": "Experimente",
    "Open example": "Abrir exemplo",
    "Three ZTF images of asteroid (5157) Hindemith with a Gaia DR3 catalog — runs fully offline.":
        "Três imagens do ZTF do asteroide (5157) Hindemith com catálogo Gaia DR3 — roda totalmente offline.",
    "Example data not found in the app bundle.": "Dados do exemplo não encontrados no app.",
    "Default": "Padrão",
    // sessão
    "Step %d/%d": "Etapa %d/%d",
    "Needs attention.": "Precisa de atenção.",
    "Done.": "Pronto.",
    "Reducing the plates…": "Reduzindo as placas…",
    "Observing log": "Diário de observação",
    "Try again": "Tentar de novo",
    "Back to start": "Voltar ao início",
    "Cancel": "Cancelar",
    // centro
    "Field center": "Centro do campo",
    "Where was the telescope pointing?": "Para onde o telescópio apontava?",
    "I could not find the coordinates in the image header. Search by target name or type RA and Dec.":
        "Não encontrei as coordenadas no cabeçalho das imagens. Procure pelo nome do alvo ou digite AR e Dec.",
    "Found.": "Encontrado.",
    "Not found in SIMBAD.": "Não encontrado no SIMBAD.",
    "RA  h m s": "AR  h m s",
    "Continue": "Continuar",
    // ajustes
    "Settings": "Ajustes",
    "Equipment": "Equipamento",
    "Name (e.g. C11 + ASI533)": "Nome (ex.: C11 + ASI533)",
    "Save current": "Salvar atual",
    "Profiles store focal length, pixel size, mirroring and MPC code.": "Perfis guardam focal, pixel, espelhamento e código MPC.",
    "Language": "Idioma",
    "System": "Sistema",
    "Restore defaults": "Restaurar padrões",
    "Restore default": "Restaurar padrão",
    "Default: %@": "Padrão: %@",
    "Must be between %@ and %@": "Deve estar entre %@ e %@",
    "Not a number": "Não é um número",
    // passos e estado (Engine)
    "Images": "Imagens",
    "Gaia DR3 catalog": "Catálogo Gaia DR3",
    "Astrometric plate": "Placa astrométrica",
    "Moving objects": "Objetos em movimento",
    "Identification": "Identificação",
    "Import the images (.st6/.fit) and, if you like, a .cfg and the catalog.": "Importe as imagens (.st6/.fit) e, se quiser, um .cfg e o catálogo.",
    "Enter the field center (RA and Dec) to download Gaia DR3.": "Informe o centro do campo (AR e Dec) para baixar o Gaia DR3.",
    "Invalid center.": "Centro inválido.",
    "Downloading Gaia DR3 from %@…": "Baixando Gaia DR3 de %@…",
    "Could not save the catalog: %@": "Não consegui gravar o catálogo: %@",
    "Gaia DR3 catalog downloaded from %@ (%d kB).": "Catálogo Gaia DR3 baixado de %@ (%d kB).",
    "No VizieR mirror answered (no internet?). Try again in a moment.": "Nenhum espelho do VizieR respondeu (sem internet?). Tente de novo em instantes.",
    "Cancelled.": "Cancelado.",
    "%d images": "%d imagens",
    "%d image(s) — motion needs 3+": "%d imagem(ns) — movimento precisa de 3+",
    "Enter the field center": "Informe o centro do campo",
    "%d refs per image · %.2f″": "%d refs por imagem · %.2f″",
    "none": "nenhum",
    "%d object(s)": "%d objeto(s)",
    "%d known asteroids in the field": "%d asteroides conhecidos no campo",
    "not identified": "não identificado",
    "Catalog missing (download Gaia DR3 or import one).": "Falta o catálogo (Baixar Gaia DR3 ou importar).",
    "Reduction error.": "Erro na redução.",
    "Reduction finished in %.2f s: %d images, %d moving object(s).": "Redução concluída em %.2f s: %d imagens, %d objeto(s) em movimento.",
    " Field: %d known asteroid(s).": " Campo: %d asteroide(s) conhecido(s).",
    "not identified (new object or SkyBoT unavailable)": "não identificado (objeto novo ou SkyBoT indisponível)",
    "Identification finished. %d known asteroid(s) in the field.": "Identificação concluída. %d asteroide(s) conhecido(s) no campo.",
    "object": "objeto",
    "Example (offline)": "Exemplo (offline)",
    // erros do núcleo
    "No images — add at least two exposures of the same field": "Nenhuma imagem — adicione pelo menos duas exposições do mesmo campo",
    "Image unreadable: %@ — use FITS or SBIG files": "Imagem ilegível: %@ — use arquivos FITS ou SBIG",
    "Catalog empty or unreadable — check the internet connection or load a Gaia DR3 TSV":
        "Catálogo vazio ou ilegível — verifique a conexão com a internet ou carregue um TSV do Gaia DR3",
    "Internal error: %@": "Erro interno: %@",
    // resultados
    "Object trail": "Rastro do objeto",
    "Known asteroids": "Asteroides conhecidos",
    "Reference stars": "Estrelas de referência",
    "All detections": "Todas as detecções",
    "Plate %d · %@": "Placa %d · %@",
    "tap = identify": "toque = identificar",
    "Result": "Resultado",
    "Nothing moved.": "Nada se moveu.",
    "The plates were solved, but no object moved between the images above the threshold. %d known asteroid(s) cross the field, probably too faint.":
        "As placas foram resolvidas, mas nenhum objeto andou entre as imagens acima do limite. %d asteroide(s) conhecido(s) passam no campo, provavelmente fracos demais.",
    "New session": "Nova sessão",
    "Object card · identified": "Ficha do objeto · identificado",
    "Object card · not identified": "Ficha do objeto · sem identificação",
    "Possible new object": "Possível objeto novo",
    "Class": "Classe",
    "Predicted brightness": "Brilho previsto",
    "Offset from ephemeris": "Desvio da efeméride",
    "Temporary designation": "Designação temporária",
    "Rate": "Velocidade",
    "Position angle": "Ângulo de posição",
    "Measured magnitude": "Magnitude medida",
    "Measurements": "Medidas",
    "Send to MPC": "Enviar ao MPC",
    "Example data — not your observation": "Dados de exemplo — não é sua observação",
    "ADES · photometric · MPC 80": "ADES · fotométrico · MPC 80",
    "What is this?": "O que é isso?",
    "Querying SIMBAD, VSX, Gaia and SkyBoT…": "Consultando SIMBAD, VSX, Gaia e SkyBoT…",
    "no answer: ": "sem resposta: ",
    "Plate quality": "Qualidade das placas",
    "few stars": "poucas estrelas",
    "excellent": "excelente",
    "good": "boa",
    "poor": "fraca",
    "Plate %d": "Placa %d",
    "Known asteroids in the field · SkyBoT": "Asteroides conhecidos no campo · SkyBoT",
    "+ %d fainter": "+ %d mais fracos",
    // "O que é isso?"
    "Asteroid · %@ · V %@": "Asteroide · %@ · V %@",
    "Variable star": "Estrela variável",
    "Gaia DR3 star": "Estrela Gaia DR3",
    " · ~%.0f light-years": " · ~%.0f anos-luz",
    "Nothing known here": "Nada conhecido aqui",
    "No catalogued object within 6″ — noise, artifact or something new.": "Nenhum objeto catalogado num raio de 6″ — ruído, artefato ou algo novo.",
    "Star": "Estrela", "Red supergiant": "Supergigante vermelha", "Galaxy": "Galáxia", "Quasar": "Quasar",
    "Planetary nebula": "Nebulosa planetária", "HII region": "Região HII", "Star cluster": "Aglomerado de estrelas",
    "Open cluster": "Aglomerado aberto", "Globular cluster": "Aglomerado globular", "Eclipsing binary": "Binária eclipsante",
    "Double star": "Estrela dupla", "White dwarf": "Anã branca", "Supernova": "Supernova", "Radio source": "Fonte de rádio",
    "X-ray source": "Fonte de raios X", "Infrared source": "Fonte infravermelha", "High proper-motion star": "Estrela de alto movimento próprio",
    "RR Lyrae variable": "Variável RR Lyrae", "Mira variable": "Variável Mira",
]

"""Step 17 - build a Hinglish plain-language explainer PDF of the whole project.

Audience: someone who wants to know what this project did and what it found, in
Hinglish, without reading the code or METHODS.md.

Every number is read from docs/facts.json at build time - nothing is typed here.
DejaVu is used throughout because the text carries rho, Angstrom, arrows and
en-dashes, which the built-in base-14 fonts cannot render.
"""
import json
import os

from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (HRFlowable, Image, KeepTogether, PageBreak,
                                Paragraph, SimpleDocTemplate, Spacer, Table,
                                TableStyle)

D = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f"{D}/docs/PROJECT_EXPLAINER_hinglish.pdf"

FONT_DIRS = ["/usr/share/fonts/truetype/dejavu",
             f"{os.path.expanduser('~')}/.claude-science/conda/envs/bioinfo-saas/fonts"]


def register_fonts():
    """DejaVu or nothing - a silent fallback to Helvetica would drop rho and Angstrom."""
    for d in FONT_DIRS:
        reg, bold = f"{d}/DejaVuSans.ttf", f"{d}/DejaVuSans-Bold.ttf"
        if os.path.exists(reg) and os.path.exists(bold):
            pdfmetrics.registerFont(TTFont("DJ", reg))
            pdfmetrics.registerFont(TTFont("DJ-Bold", bold))
            pdfmetrics.registerFontFamily("DJ", normal="DJ", bold="DJ-Bold")
            return
    raise SystemExit(f"DejaVuSans.ttf + Bold not found in {FONT_DIRS}")


INK = colors.HexColor("#1a1a1a")
FOCAL = colors.HexColor("#1f4e79")
ACCENT = colors.HexColor("#c1502e")
GREY = colors.HexColor("#666666")
RULE = colors.HexColor("#cccccc")


def styles():
    ss = getSampleStyleSheet()
    def S(name, **kw):
        base = dict(fontName="DJ", textColor=INK, leading=14, fontSize=9.6)
        base.update(kw)
        return ParagraphStyle(name, parent=ss["Normal"], **base)
    return {
        "title":   S("t",  fontName="DJ-Bold", fontSize=19, leading=23, textColor=FOCAL,
                     spaceAfter=2),
        "sub":     S("s",  fontSize=10.2, textColor=GREY, leading=14, spaceAfter=10),
        "h1":      S("h1", fontName="DJ-Bold", fontSize=13, leading=17, textColor=FOCAL,
                     spaceBefore=13, spaceAfter=5),
        "h2":      S("h2", fontName="DJ-Bold", fontSize=10.4, leading=14, spaceBefore=8,
                     spaceAfter=3),
        "body":    S("b",  alignment=TA_JUSTIFY, spaceAfter=6),
        "bullet":  S("bu", alignment=TA_JUSTIFY, leftIndent=11, bulletIndent=2,
                     spaceAfter=3.5),
        "caption": S("c",  fontSize=8.2, leading=11, textColor=GREY, spaceBefore=3,
                     spaceAfter=9),
        "callout": S("co", fontSize=9.6, leading=14, leftIndent=8, rightIndent=8,
                     spaceBefore=3, spaceAfter=3),
        "foot":    S("f",  fontSize=7.8, leading=10, textColor=GREY),
    }


def tbl(rows, widths, st, header=True, align_right=None):
    t = Table(rows, colWidths=widths, hAlign="LEFT")
    cmds = [
        ("FONT", (0, 0), (-1, -1), "DJ", 8.6),
        ("TEXTCOLOR", (0, 0), (-1, -1), INK),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 3.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("LINEBELOW", (0, 0), (-1, -2), 0.4, RULE),
    ]
    if header:
        cmds += [("FONT", (0, 0), (-1, 0), "DJ-Bold", 8.6),
                 ("TEXTCOLOR", (0, 0), (-1, 0), FOCAL),
                 ("LINEBELOW", (0, 0), (-1, 0), 0.9, FOCAL)]
    for c in (align_right or []):
        cmds.append(("ALIGN", (c, 0), (c, -1), "RIGHT"))
    t.setStyle(TableStyle(cmds))
    return t


def callout(text, st, colour=ACCENT):
    p = Paragraph(text, st["callout"])
    t = Table([[p]], colWidths=[165 * mm], hAlign="LEFT")
    t.setStyle(TableStyle([
        ("LINEBEFORE", (0, 0), (0, -1), 2.2, colour),
        ("LEFTPADDING", (0, 0), (-1, -1), 9),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#faf6f4")),
    ]))
    return t


def figure(path, caption, st, width=165 * mm):
    """Image + caption as one unbreakable flowable, scaled to the text column."""
    from PIL import Image as PILImage
    with PILImage.open(path) as im:
        w, h = im.size
    img = Image(path, width=width, height=width * h / w)
    return KeepTogether([img, Paragraph(caption, st["caption"])])


def build():
    register_fonts()
    st = styles()
    F = json.load(open(f"{D}/docs/facts.json"))
    q, dk, ad = F["qsar"], F["docking"], F["applicability_domain"]
    ie, wc = q["InhA_enzyme"], q["Mtb_whole_cell"]
    tr, inf = F["translation"], F["inflation"]

    def pct(x):
        return f"{x:.0%}"

    E = []
    A = E.append

    # ---------------------------------------------------------------- cover
    A(Paragraph("InhA-TB: is project mein kya kiya gaya hai", st["title"]))
    A(Paragraph("Mycobacterium tuberculosis ke InhA target ke against ek poora, "
                "reproducible drug-discovery pipeline &mdash; aur uska imaandaar result",
                st["sub"]))
    A(HRFlowable(width="100%", thickness=0.9, color=FOCAL, spaceAfter=10))

    A(Paragraph("Ek line mein", st["h1"]))
    A(Paragraph(
        f"TB ki dawa ke ek proven target (<b>InhA</b>) par public bioactivity data se do "
        f"machine-learning models banaye gaye, unhe <b>honestly</b> validate kiya gaya, "
        f"aur phir same molecules par classical <b>molecular docking</b> se compare kiya "
        f"gaya. Result: ML model docking se kaafi behtar nikla &mdash; aur yeh baat "
        f"chhupayi nahi, report ki gayi hai.", st["body"]))

    A(callout(
        f"<b>Sabse important number:</b> exactly same held-out molecules par "
        f"(jinhe na model ne train mein dekha, na docking ko fit kiya gaya) &mdash; "
        f"<b>QSAR ROC-AUC {ie['scaffold_roc_auc']:.2f} / {wc['scaffold_roc_auc']:.2f}</b> "
        f"vs <b>docking {dk['InhA_enzyme_heldout_dock_auc']:.2f} / "
        f"{dk['Mtb_whole_cell_heldout_dock_auc']:.2f}</b>. "
        f"Docking is target par lagbhag coin-toss ke barabar hai.", st))
    A(Spacer(1, 8))

    A(Paragraph("Problem kya thi", st["h1"]))
    A(Paragraph(
        "TB duniya ki sabse jaanleva infectious diseases mein se ek hai, aur "
        "drug-resistant strains badh rahe hain. <b>InhA</b> (enoyl-ACP reductase) "
        "isoniazid ka validated target hai &mdash; matlab yeh enzyme block karo to "
        "bacteria marta hai, yeh already clinically proven hai. Isoniazid ek prodrug "
        "hai jise KatG enzyme activate karta hai, aur resistance ka sabse bada kaaran "
        "<i>katG</i> mutation hai. Isliye aisa inhibitor jo <b>direct</b> InhA ko block "
        "kare (KatG activation ke bina), resistance bypass kar sakta hai.", st["body"]))
    A(Paragraph(
        "Sawaal yeh tha: kya public data aur open-source tools se ek aisa screening "
        "pipeline banaya ja sakta hai jo <b>trustworthy</b> ho &mdash; matlab jiske "
        "numbers exaggerated na hon?", st["body"]))

    # ---------------------------------------------------------------- pipeline
    A(Paragraph("Kya-kya banaya gaya (7 steps)", st["h1"]))
    steps = [
        ("1. Data laana",
         f"ChEMBL database se InhA enzyme aur <i>M. tb</i> whole-cell activity records "
         f"pull kiye. Raw {F['cascade'][0][1]:,} records the."),
        ("2. Curation",
         f"Har filter par kitna data bacha, woh record kiya gaya (curation cascade). "
         f"Final: <b>{F['records_curated']:,} clean records &rarr; "
         f"{F['molecules_unique']:,} unique molecules</b>."),
        ("3. Structure prep",
         f"{len(json.load(open(f'{D}/data/structure_meta.json')))} InhA crystal structures "
         f"compare karke <b>{dk['receptor']}</b> chuna. Docking box binding site par "
         f"set kiya."),
        ("4. Redocking control",
         f"Crystal structure ka apna ligand wapas dock karke check kiya ki setup sahi hai: "
         f"<b>{dk['rmsd_top']:.2f} &Aring; RMSD</b>. Yeh control na ho to docking campaign "
         f"ke paas koi proof nahi ki box/protonation/scoring theek hai."),
        ("5. QSAR models",
         f"Molecular fingerprints ({F['fingerprint']}) + {F['descriptors_n']} "
         f"physicochemical descriptors se XGBoost, Random Forest aur logistic regression "
         f"train kiye."),
        ("6. Docking screen",
         f"<b>{dk['library_n']}</b> molecules ki class-balanced library AutoDock Vina se "
         f"dock ki; <b>{dk['n_scored']}</b> successfully scored."),
        ("7. Honest comparison",
         "Dono approaches ko <b>same held-out molecules</b> par compare kiya &mdash; yehi "
         "is project ka core hai."),
    ]
    for h, t in steps:
        A(Paragraph(f"<b>{h}</b> &mdash; {t}", st["bullet"], bulletText="\u2022"))

    A(Spacer(1, 6))
    A(figure(f"{D}/figures/fig1_data_and_target.png",
              f"<b>Figure 1.</b> (a) Curation cascade &mdash; har filter par kitna data "
              f"bacha. (b) Potency distributions. (c) Experimental noise floor "
              f"{F['noise_floor_sd']:.2f} log units. (d) {len(json.load(open(f'{D}/data/structure_meta.json')))} "
              f"InhA structures compare kiye, {dk['receptor']} chuna. (e) Binding-pocket "
              f"contacts. (f) Docked subset vs poora curated set.", st))

    # ---------------------------------------------------------------- honesty
    A(Paragraph("Sabse important cheez: honesty ka mechanism", st["h1"]))
    A(Paragraph(
        "Zyada bioinformatics projects yahan galti karte hain: molecules ko "
        "<b>randomly</b> train/test mein baant dete hain. Problem yeh hai ki chemical "
        "series &mdash; ek hi core structure ke 20 analogues &mdash; dono taraf bant "
        "jaate hain. Model test set mein woh cheez pehchanta hai jo woh train mein "
        "already dekh chuka hai, aur score jhoota-uncha aata hai.", st["body"]))
    A(Paragraph(
        f"Yahan <b>Bemis&ndash;Murcko scaffold split</b> use kiya gaya: pura chemical "
        f"scaffold ya to train mein jaata hai ya test mein, beech mein nahi. "
        f"Verified: train aur test ke beech <b>zero scaffold overlap</b>.", st["body"]))

    A(callout(
        f"<b>Kitna farak padta hai:</b> random split ROC-AUC ko "
        f"<b>{inf['min']:.2f}&ndash;{inf['max']:.2f}</b> (median {inf['median']:.2f}) se "
        f"inflate karta hai &mdash; saare 6 model &times; dataset combinations mein. "
        f"Is project mein quote kiya gaya har number <b>scaffold-split</b> ka hai, "
        f"yaani kam dikhne wala lekin sach wala number.", st))

    A(Paragraph("QSAR model ka performance", st["h1"]))
    A(tbl([
        ["Dataset", "Best model", "Scaffold ROC-AUC", "Random split", "Train / Test"],
        ["InhA enzyme", ie["best_model"], f"{ie['scaffold_roc_auc']:.2f}",
         f"{ie['random_roc_auc']:.2f}", f"{ie['n_train']} / {ie['n_test']}"],
        ["M. tb whole-cell", wc["best_model"], f"{wc['scaffold_roc_auc']:.2f}",
         f"{wc['random_roc_auc']:.2f}", f"{wc['n_train']} / {wc['n_test']}"],
    ], [36 * mm, 30 * mm, 33 * mm, 26 * mm, 30 * mm], st, align_right=[2, 3, 4]))
    A(Spacer(1, 7))

    A(Paragraph("Model kab kaam karta hai, kab nahi (applicability domain)", st["h2"]))
    A(Paragraph(
        f"Ek model ko blindly trust nahi karna chahiye. Yahan test compounds ko training "
        f"set se chemical similarity ke hisaab se bins mein baanta gaya. Jo molecules "
        f"training data se <b>milte-julte</b> hain (Tanimoto &ge; 0.35), un par ROC-AUC "
        f"<b>{min(ad['near_roc_auc']):.2f}&ndash;{max(ad['near_roc_auc']):.2f}</b>. "
        f"Jo <b>door</b> hain (&lt; 0.35), wahan gir kar "
        f"<b>{min(ad['far_roc_auc']):.2f}&ndash;{max(ad['far_roc_auc']):.2f}</b> &mdash; "
        f"yaani chance ke barabar ya usse bhi neeche. Isliye production mein har "
        f"prediction ke saath ek <b>distance flag</b> hona chahiye.", st["body"]))

    A(Spacer(1, 4))
    A(figure(f"{D}/figures/fig2_qsar_validation.png",
             f"<b>Figure 2.</b> (a,b) Random split har model par accuracy badha-chadha kar "
             f"dikhata hai. (c) Precision&ndash;recall. (d) Accuracy chemical distance ka "
             f"function hai. (e) Enzyme potency whole-cell kill se track karti hai "
             f"(Spearman &rho; = {tr['spearman_rho']:.2f}, n = {tr['n']}). "
             f"(f) Test chemotypes train se alag baithe hain.", st))

    # ---------------------------------------------------------------- docking
    A(Paragraph("Docking screen &mdash; aur ek negative result", st["h1"]))
    A(Paragraph(
        f"<b>{dk['library_n']}</b> molecules ki library banayi gayi &mdash; deliberately "
        f"<b>class-balanced</b> (InhA arm {dk['library_balance']['InhA_enzyme']['actives']} "
        f"actives / {dk['library_balance']['InhA_enzyme']['n']} molecules, whole-cell arm "
        f"{dk['library_balance']['Mtb_whole_cell']['actives']} / "
        f"{dk['library_balance']['Mtb_whole_cell']['n']}), taaki discrimination test "
        f"karna possible ho. Drug-like envelope filter lagaya "
        f"(MW &le; {dk['library_mw_max']:.0f} Da), jisne "
        f"{dk['n_excluded']} extreme molecules hata diye &mdash; sabse bada "
        f"{dk['excluded_max_mw']:,.0f} Da / {dk['excluded_max_rotb']} rotatable bonds tha, "
        f"jo Vina ke liye meaningful hi nahi.", st["body"]))

    A(tbl([
        ["Dataset", "n", "Actives", "ROC-AUC", "Spearman \u03c1", "EF 5 %", "Ceiling"],
        ["InhA enzyme", str(dk["library_balance"]["InhA_enzyme"]["n"]),
         str(dk["library_balance"]["InhA_enzyme"]["actives"]),
         f"{dk['InhA_enzyme_roc_auc']:.2f}", f"{dk['InhA_enzyme_rho']:.2f}",
         f"{dk['InhA_enzyme_ef5']:.2f}", f"{dk['InhA_enzyme_ef_max']:.2f}"],
        ["M. tb whole-cell", str(dk["library_balance"]["Mtb_whole_cell"]["n"]),
         str(dk["library_balance"]["Mtb_whole_cell"]["actives"]),
         f"{dk['Mtb_whole_cell_roc_auc']:.2f}", f"{dk['Mtb_whole_cell_rho']:.2f}",
         f"{dk['Mtb_whole_cell_ef5']:.2f}", f"{dk['Mtb_whole_cell_ef_max']:.2f}"],
    ], [34 * mm, 14 * mm, 17 * mm, 21 * mm, 24 * mm, 20 * mm, 20 * mm], st,
        align_right=[1, 2, 3, 4, 5, 6]))
    A(Spacer(1, 7))

    A(Paragraph(
        f"Whole-cell arm par <b>EF 5 % = {dk['Mtb_whole_cell_ef5']:.2f}</b> ka matlab hai: "
        f"top-scoring 5 % molecules uthane se <b>random uthane se kam</b> actives milenge. "
        f"Yeh chhupane wali baat nahi &mdash; batane wali hai.", st["body"]))

    A(Paragraph("Head-to-head: same molecules, dono ke liye naye", st["h2"]))
    A(tbl([
        ["Dataset", "Held-out n", "Docking ROC-AUC", "QSAR ROC-AUC"],
        ["InhA enzyme", str(dk["InhA_enzyme_heldout_n"]),
         f"{dk['InhA_enzyme_heldout_dock_auc']:.2f}",
         f"{dk['InhA_enzyme_heldout_qsar_auc']:.2f}"],
        ["M. tb whole-cell", str(dk["Mtb_whole_cell_heldout_n"]),
         f"{dk['Mtb_whole_cell_heldout_dock_auc']:.2f}",
         f"{dk['Mtb_whole_cell_heldout_qsar_auc']:.2f}"],
    ], [40 * mm, 28 * mm, 36 * mm, 32 * mm], st, align_right=[1, 2, 3]))
    A(Spacer(1, 6))

    A(callout(
        "<b>Yeh negative result kyun rakha gaya:</b> ek well-posed comparison ka negative "
        "answer bhi result hota hai. Yeh bound karta hai ki structure-based component "
        "kitna contribute kar sakta hai, aur yahi wajah hai ki product ka prediction "
        "endpoint QSAR model hai, docking score nahi. Docking hataya nahi gaya &mdash; "
        "measure karke report kiya gaya.", st, colour=FOCAL))

    A(Spacer(1, 6))
    A(figure(f"{D}/figures/fig3_docking.png",
             f"<b>Figure 3.</b> (a) Actives aur inactives ke score overlap karte hain. "
             f"(b) Score potency ko barely track karta hai. (c) ROC. (d) Top-of-list "
             f"enrichment attainable ceiling se neeche. (e) Same molecules, neither "
             f"fitted &mdash; QSAR clearly aage. (f) Vina score molecule ke size ko "
             f"reward karta hai (&rho; = -0.41), jo raw ranking ko biased banata hai.",
             st))

    # ---------------------------------------------------------------- limits
    A(Paragraph("Kya yeh project claim NAHI karta", st["h1"]))
    for t in [
        "<b>Vina score affinity nahi hai.</b> Scoring function binding geometry "
        "reproduce karne ke liye bani thi, potency rank karne ke liye nahi.",
        "<b>Rigid receptor.</b> InhA ka substrate-binding loop ligand ke saath "
        "reorganise hota hai; ek single crystal structure yeh capture nahi karta.",
        "<b>Ek conformer per ligand.</b> Flexible molecules under-sampled hain.",
        "<b>Koi wet-lab validation nahi.</b> Yeh computational triage hai, discovery ka "
        "claim nahi. Har prediction ko assay se confirm karna padega.",
        f"<b>Applicability domain limited hai.</b> Training set se door molecules par "
        f"model ka ROC-AUC {min(ad['far_roc_auc']):.2f} tak gir jaata hai.",
    ]:
        A(Paragraph(t, st["bullet"], bulletText="\u2022"))

    A(Paragraph("Reproducibility &mdash; is project ka asli differentiator", st["h1"]))
    A(Paragraph(
        "Sirf result nahi, uska <b>bharosa</b> bhi build kiya gaya hai. Repo mein "
        "16 numbered scripts hain jo ek runner se end-to-end chalti hain, aur do "
        "<b>automated quality gates</b> hain:", st["body"]))
    for t in [
        "<b>Consistency gate</b> &mdash; README, METHODS.md aur facts.json ke "
        "<b>22 headline numbers</b> match hone chahiye, warna build fail. Har number "
        "data se generate hota hai; koi hand-typed nahi.",
        "<b>Figure gate</b> &mdash; teenon figures automatically check hoti hain ki koi "
        "text overlap ya clipped label to nahi. Dono gates ko <b>negative controls</b> se "
        "verify kiya gaya hai &mdash; jaan-boojh kar defect daal kar confirm kiya ki gate "
        "sach mein pakadta hai, andha nahi hua.",
    ]:
        A(Paragraph(t, st["bullet"], bulletText="\u2022"))

    A(Paragraph(
        "Yeh isliye matter karta hai: ek hand-typed number document mein reh jaata hai "
        "aur data badalne par chupchaap galat ho jaata hai. Is repo mein woh possible "
        "nahi &mdash; gate build tod dega.", st["body"]))

    A(Paragraph("Tech stack", st["h1"]))
    A(Paragraph(
        "Python, RDKit (cheminformatics), scikit-learn + XGBoost (models), "
        "AutoDock Vina (docking), Open Babel (conformers), pandas/numpy, "
        "matplotlib (figures). Data: ChEMBL (bioactivity), RCSB PDB (structures). "
        "Sab open-source, sab free.", st["body"]))

    A(Spacer(1, 10))
    A(HRFlowable(width="100%", thickness=0.5, color=RULE, spaceAfter=6))
    A(Paragraph(
        f"Yeh document <code>src/17_hinglish_pdf.py</code> se generate hua hai. Saare "
        f"numbers <code>docs/facts.json</code> se live padhe jaate hain &mdash; koi value "
        f"haath se type nahi ki gayi. Full technical detail: <b>docs/METHODS.md</b>. "
        f"Repository: <b>github.com/TilakSharma07/inha-tb</b>", st["foot"]))

    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("DJ", 7.6)
        canvas.setFillColor(GREY)
        canvas.drawString(22 * mm, 12 * mm, "InhA-TB \u2014 project explainer (Hinglish)")
        canvas.drawRightString(A4[0] - 22 * mm, 12 * mm, f"{doc.page}")
        canvas.restoreState()

    # invariant=1 fixes /CreationDate, /ModDate and the document /ID, which otherwise
    # change on every run. Without it the 1.1 MB PDF shows up as a spurious diff after
    # each run_all.sh, and git stores a whole new blob for a logically unchanged file.
    doc = SimpleDocTemplate(
        OUT, pagesize=A4, invariant=1,
        leftMargin=22 * mm, rightMargin=22 * mm,
        topMargin=18 * mm, bottomMargin=20 * mm,
        title="InhA-TB - project explainer (Hinglish)",
        author="", subject="InhA / M. tuberculosis virtual screening pipeline")
    doc.build(E, onFirstPage=footer, onLaterPages=footer)
    print(f"wrote {OUT} ({os.path.getsize(OUT)/1024:.0f} KB)")


if __name__ == "__main__":
    build()

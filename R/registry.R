# registry.R -- the tool registry, and the two things every tool's right rail
# is built from.
#
# This exists so app.R can stay a shell. TOOLS is the single list that drives
# the left nav, the tabsetPanel, the "Quick actions" links, AND (since the
# decomposition) which R/panels/<id>.R files get sourced -- so adding a tool is
# one entry here plus one file, with no third place to forget.
#
# The invariant that matters: every id in TOOLS must have a matching
# R/panels/<id>.R defining panel_<id>() and server_<id>(). app.R's source loop
# fails loudly if it doesn't, which is the behaviour you want -- a tool that
# silently doesn't load is much worse than one that stops the app at startup.
#
# TOOL_ABOUT / TOOL_RELATED are kept here rather than next to each panel on
# purpose: they are cross-tool editorial copy (a one-line honest description,
# and the 2 tools a typical workflow chains into), and they read better -- and
# stay consistent -- when the whole set is visible in one place.

`%||%` <- function(x, default) if (is.null(x) || (length(x) == 1 && is.na(x)) || (is.character(x) && !nzchar(x))) default else x

TOOLS <- list(
  list(id = "home",     label = "Home",              icon = "\U0001f3e0", group = "Home"),
  list(id = "notebook", label = "Lab Notebook",      icon = "\U0001f4d3", group = "Notebook"),
  list(id = "explorer", label = "Transcript Explorer", icon = "\U0001f9ed", group = "Design"),
  list(id = "extractor", label = "Exon Extractor",     icon = "\U00002702", group = "Design"),
  list(id = "design",   label = "Primer & Schematic", icon = "\U0001f9ec", group = "Design"),
  list(id = "cryptic",  label = "Cryptic Splicing Engine", icon = "\U0001f52c", group = "Design"),
  list(id = "splicecode", label = "Splice Code",        icon = "\U0001f9ee", group = "Design"),
  list(id = "consequence", label = "Protein Consequence", icon = "\U0001f9e9", group = "Design"),
  list(id = "batch",    label = "Panel Runner",       icon = "\U0001f5c2", group = "Design"),
  list(id = "plasmid",  label = "Plasmid Creator",    icon = "\U0001f504", group = "Design"),
  list(id = "gibson",   label = "Gibson Assembly",    icon = "\U0001f517", group = "Design"),
  list(id = "plasmidqc", label = "Plasmid QC",         icon = "\U0001f50e", group = "Design"),
  list(id = "pcr",      label = "PCR Setup",          icon = "\U0001f9ea", group = "Design"),
  list(id = "qpcr",     label = "qPCR (ddCt)",        icon = "\U0001f4c9", group = "Analysis"),
  list(id = "dens",     label = "Densitometry",       icon = "\U0001f4ca", group = "Analysis"),
  list(id = "sc",       label = "Standard Curve",     icon = "\U0001f4c8", group = "Analysis"),
  list(id = "report",   label = "Methods & Ordering", icon = "\U0001f4cb", group = "Analysis"),
  list(id = "norm",     label = "Protein Normalization", icon = "⚖️", group = "Calculators"),
  list(id = "a280",     label = "A280 Calculator",    icon = "\U0001f9eb", group = "Calculators"),
  list(id = "pp",       label = "Protein Parameters", icon = "\U0001f9ec", group = "Calculators"),
  list(id = "dil",      label = "Dilution Calculator", icon = "\U0001f4a7", group = "Calculators"),
  list(id = "about",    label = "About",              icon = "\U0001f464", group = "About")
)
TOOL_BY_ID <- setNames(TOOLS, sapply(TOOLS, `[[`, "id"))

# Landing-page feature cards -- each jumps into the named tool (see server).
HOME_FEATURES <- list(
  list(tool = "cryptic",  icon = "\U0001f52c", title = "Cryptic Splicing Engine",
       blurb = "Read control vs. knockdown BAMs over a locus and flag cryptic exons, cryptic splice sites, exitrons, and intron retention in an IGV-style sashimi plot."),
  list(tool = "design",   icon = "\U0001f9ec", title = "Primer & Schematic",
       blurb = "Design junction-spanning primers against live reference sequence, with expected canonical vs. cryptic product sizes and a schematic."),
  list(tool = "notebook", icon = "\U0001f4d3", title = "Lab Notebook",
       blurb = "Build reusable procedures and spin up experiments with editable tables — everything saves to disk and reopens across sessions."),
  list(tool = "batch",    icon = "\U0001f5c2", title = "Panel Runner",
       blurb = "Run the cryptic-detection pipeline across a whole gene list against one BAM pair, then click any hit to open it in the engine."),
  list(tool = "qpcr",     icon = "\U0001f4c9", title = "qPCR & Densitometry",
       blurb = "Relative expression by 2^-ΔΔCt with QC warnings, plus Western-band densitometry normalized to a loading control."),
  list(tool = "report",   icon = "\U0001f4cb", title = "Methods & Ordering",
       blurb = "Auto-assemble a primer ordering sheet, a references list, and a templated Methods paragraph from whatever you ran this session.")
)

# right-rail content: a one-line honest description (no invented metrics) and
# 2 related tools to jump to -- both genuinely derived from how the tools chain
# together in a typical workflow.
TOOL_ABOUT <- list(
  home     = "The landing page: what Lit2Bench does and who built it, with quick jumps into the main tools.",
  notebook = "An electronic lab notebook that saves to disk. Procedures are reusable templates (Objective / Reagents / Experimental setup with editable tables / Results / Conclusions); spin up an Experiment from a procedure, fill in your run and results, and Save — entries persist as JSON files you can reopen and edit across sessions.",
  explorer = "Search a gene symbol, RefSeq/Ensembl transcript ID, or a locus and see every annotated isoform in the region -- strand, length, coding status, exon/intron counts, CDS and UTR spans.",
  extractor = "Pulls real exon and intron sequence for a chosen transcript, exports BED/FASTA/CSV/JSON/GTF, and hands a chosen exon straight to the Primer Designer -- no manual coordinate copying.",
  design  = "Designs a junction-spanning primer pair against live reference sequence and shows the expected canonical vs. cryptic-exon product sizes.",
  cryptic = "Reads control vs. knockdown BAMs over a locus and screens for all four recognized types of splicing defect -- cryptic exon inclusion, cryptic splice site selection, exitrons, and intron retention -- via an IGV-style sashimi plot without leaving the app. Candidate exons carry coverage evidence as well as junction evidence: whether the exon body itself is elevated in knockdown, scored against its own intron. Use it as a screen for which loci deserve a look -- the confidence tier grades the flanking junctions, not the span, so check coordinates on the plot before designing against them.",
  splicecode = "Answers why an exon would be silent in the first place: scores every splice site in a transcript against a matrix built from real annotated sites, so you can see where a cryptic site ranks among the gene's own; then reads the polypyrimidine tract, the branch-point consensus, and two independent kinds of RNA-binding-protein evidence -- strand-aware motif richness, and measured binding from published ENCODE eCLIP. Pick from 176 proteins (TDP-43, RBFOX2, PTBP1, hnRNP C, NOVA, MBNL1, QKI and more); only the binding evidence depends on the choice. Needs no BAMs; a gene symbol is enough.",
  consequence = "Takes a cryptic exon and works out what it does to the protein: splices it into the real annotated transcript, translates from the annotated start codon, and reports the frameshift, the premature stop, whether the transcript is predicted to be degraded by NMD (55-nt rule), and which UniProt domains the truncation removes.",
  batch   = "Runs the Cryptic Splicing Engine's full detection across a whole list of genes/loci against one BAM pair, returning a row-per-locus summary. Click any row to open that locus in the engine.",
  plasmid = "Joins your parts end-to-end, circularizes them, and draws the resulting map.",
  gibson  = "Designs the primer pairs to Gibson/NEBuilder-assemble an ordered set of DNA fragments — a gene-specific annealing region sized to your target Tm plus the homology tail that builds each junction, for a circular (insert+vector) or linear assembly.",
  plasmidqc = "QC your Plasmidsaurus (or any) sequencing reads against a reference plasmid: local pairwise alignment on both strands, a PASS / GENE_FOUND / FLAGGED verdict from identity + coverage thresholds, called substitutions/indels and truncations, and unmatched flanks you can screen against NCBI. Ported from Alex Luu's GeneAlign.",
  pcr     = "Scales stock/final concentrations into a master mix for N reactions plus excess.",
  qpcr    = "Calculates ΔΔCt relative expression against a chosen calibrator sample.",
  dens    = "Normalizes target band intensity to a loading control across lanes.",
  sc      = "Fits a standard curve (linear or quadratic) and back-calculates unknown concentrations.",
  report  = "Collects everything this session produced into three paste-ready artifacts: a primer ordering sheet, a templated Methods paragraph, and a Methods references list — all built from what you actually ran, with blanks marked rather than invented.",
  norm    = "Works out lysate / water / dye volumes for equal-protein-mass loading.",
  a280    = "Converts A280 absorbance to concentration via Beer-Lambert, using your extinction coefficient and MW.",
  pp      = "Computes MW, pI, and extinction coefficient directly from an amino-acid sequence.",
  dil     = "Solves C1V1 = C2V2 for stock volume and diluent volume.",
  about   = "Who built Lit2Bench, and the lab behind the cryptic-splicing biology it targets."
)
TOOL_RELATED <- list(
  home = c("cryptic", "notebook"),
  notebook = c("design", "pcr"),
  explorer = c("extractor", "design"), extractor = c("design", "explorer"),
  design = c("extractor", "cryptic", "plasmid"), cryptic = c("splicecode", "consequence"),
  splicecode = c("cryptic", "consequence"),
  consequence = c("cryptic", "splicecode"),
  batch = c("cryptic", "design"),
  plasmid = c("gibson", "plasmidqc"), gibson = c("plasmid", "pcr"),
  plasmidqc = c("plasmid", "gibson"),
  pcr = c("design", "qpcr"),
  qpcr = c("dens", "sc"), dens = c("norm", "sc"), sc = c("norm", "a280"),
  report = c("design", "qpcr"),
  norm = c("a280", "dil"), a280 = c("pp", "dil"), pp = c("a280", "dil"), dil = c("pcr", "norm"),
  about = c("notebook", "cryptic")
)

l2b_generic_aside <- function(id, status_ui) {
  tagList(
    l2b_aside_card("About this tool",
      p(style = "font-size:13.5px; color:var(--l2b-text-muted); line-height:1.5; margin:0;", TOOL_ABOUT[[id]])),
    l2b_aside_card("Status", status_ui),
    l2b_aside_card("Quick actions",
      lapply(TOOL_RELATED[[id]], function(rid) {
        t <- TOOL_BY_ID[[rid]]
        l2b_aside_link(paste0("aside_nav_", rid), t$icon, t$label)
      }))
  )
}

# The standard right-rail status line: an error if the tool errored, a neutral
# note if it hasn't run, and a tool-written one-liner about the real result
# otherwise. Lives here (not in a panel) because every tool file calls it from
# the aside function it publishes -- see l2b_new_ctx()'s $publish in R/ctx.R.
status_row <- function(res, err, ready_fn) {
  if (!is.null(err)) return(l2b_aside_status(FALSE, err))
  if (is.null(res)) return(div(class = "l2b-aside-note", "No results yet."))
  l2b_aside_status(TRUE, ready_fn(res))
}

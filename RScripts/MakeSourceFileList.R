# SCRIPT TO GENERATE WORD FILE OF SOURCE REFERENCES

# Read in list of source files:
SourceFiles <- gsub(".xml", "", list.files(path = "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/InputData/XML"))

# Read in primates HTML:
PrimatesHTML <- c(readLines("http://www.graemetlloyd.com/matrprim.html"), "<p class=\"hangingindent\">Silcox, M. T., Bloch, J. I., Boyer, D. M. and Houde, P., 2010.", "Cranial anatomy of Paleocene and Eocene <em>Labidolemur kayi</em> (Mammalia: Apatotheria), and the relationships of the Apatemyidae to other mammals.", "<em>Zoological Journal of the Linnean Society</em>, <b>160</b>, 773-825.<br><font size=\"-1\">", "<a href=\"nexus/Silcox_etal_2010a.nex\" target=\"_blank\">NEXUS</a> | <a href=\"tnt/Silcox_etal_2010a.tnt\" target=\"_blank\">TNT</a> | <a href=\"mpts/Silcox_etal_2010a.tre\" target=\"_blank\">MPT(s)</a> <a href=\"firstmpt/Silcox_etal_2010a.tre\" target=\"_blank\">(1)</a> | <a href=\"sc/Silcox_etal_2010a.tre\" target=\"_blank\">SC</a> <a href=\"schtml/Silcox_etal_2010a.html\" target=\"_blank\">(HTML)</a> | <a href=\"mrp/Silcox_etal_2010amrp.nex\" target=\"_blank\">MRP</a> | <a href=\"xml/Silcox_etal_2010a.xml\" target=\"_blank\">XML</a></font></p>")

# Get list of HTML formatted references:
HTMLFormattedRefs <- lapply(as.list(grep("<p class=\"hangingindent\">", PrimatesHTML)), function(x) PrimatesHTML[x:length(PrimatesHTML)][1:grep("</p>", PrimatesHTML[x:length(PrimatesHTML)])[1]])

# Get unique paper titles from XML:
PaperTitles <- unique(unlist(lapply(as.list(list.files(path = "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/InputData/XML", full.names=TRUE)), function(x) metatree::ReadMetatreeXML(x)$SourceTree$Source$Title$TagContents)))

# Get just the source references and format as markdown:
MarkdownSourceRefs <- paste(unlist(lapply(as.list(PaperTitles), function(x) {
  HTMLRef <- HTMLFormattedRefs[[which(unlist(lapply(HTMLFormattedRefs, function(y) length(grep(toupper(x), toupper(y), fixed = TRUE)))) == 1)]]
  while(length(grep("  ", HTMLRef)) > 0) HTMLRef <- gsub("  ", " ", HTMLRef)
  HTMLRef <- paste(HTMLRef, collapse = "")
  HTMLRef <- strsplit(HTMLRef, split = "<br>")[[1]][1]
  HTMLRef <- gsub(" <p class=\"hangingindent\">", "", HTMLRef)
  HTMLRef <- gsub(", <b></b>, .", ".", HTMLRef)
  HTMLRef <- gsub("<b>|</b>", "**", HTMLRef)
  HTMLRef <- gsub("<em>|</em>", "*", HTMLRef)
  HTMLRef
})), sep = "\n")

# Write out markdown file:
write(paste("\\setlength{\\parindent}{-0.2in}", "\\setlength{\\leftskip}{0.2in}", "\\setlength{\\parskip}{8pt}", "\\noindent", "\n# References for included source data\n", paste(unique(MarkdownSourceRefs), collapse = "\n\n"), sep = "\n"), file = "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/InputData/SourceList.md")

# Convert markdown file to MS Word:
setwd("~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/InputData")
rmarkdown::pandoc_convert("SourceList.md", output = "SourceList.docx")

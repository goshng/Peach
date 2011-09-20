INSERT INTO dbDb
    (name, description, nibPath, organism,
     defaultPos, active, orderKey, genome, scientificName,
     htmlPath, hgNearOk, hgPbOk, sourceName)
VALUES
    ("strPyg1", "Sep 19 2011", "/gbdb/strPyg1", "SpyogenesMGAS315",
     "chr1:1-10000", 1, 1, "SpyogenesMGAS315", "Streptococcus mutans",
     "/gbdb/strPyg1/html/description.html", 0, 0, "Streptococcus Recombination study version 1.0");
INSERT INTO defaultDb (genome, name) VALUES ("SpyogenesMGAS315", "strPyg1");
INSERT INTO genomeClade (genome, clade, priority) VALUES ("SpyogenesMGAS315", "streptococcus", 2);
/* INSERT INTO clade (name, label, priority) VALUES ("streptococcus", "Streptococcus", 45); */

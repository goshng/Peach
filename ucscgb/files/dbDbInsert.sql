INSERT INTO dbDb
    (name, description, nibPath, organism,
     defaultPos, active, orderKey, genome, scientificName,
     htmlPath, hgNearOk, hgPbOk, sourceName)
VALUES
    ("strMut2", "Sep 20 2011", "/gbdb/strMut2", "SmutansUA195",
     "chr1:1-10000", 1, 1, "SmutansUA195", "Streptococcus mutans",
     "/gbdb/strMut2/html/description.html", 0, 0, "S. mutans RNA-Seq data version 1.0");
INSERT INTO defaultDb (genome, name) VALUES ("SmutansUA195", "strMut2");
INSERT INTO genomeClade (genome, clade, priority) VALUES ("SmutansUA195", "streptococcus", 1);
INSERT INTO clade (name, label, priority) VALUES ("streptococcus", "Streptococcus", 45);

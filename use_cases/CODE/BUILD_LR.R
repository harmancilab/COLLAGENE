
args = commandArgs(trailingOnly=TRUE)

DATA_MODEL_DIR <- args[1];

print(sprintf("Data dir: %s", DATA_MODEL_DIR));

feats <- read.delim(sprintf('%s/feat_matrix.txt', DATA_MODEL_DIR), header=TRUE);
pheno <- read.delim(sprintf('%s/phenotypes.txt', DATA_MODEL_DIR),header=TRUE);

filt_feats <- feats[, -c(1,2)];
filt_pheno <- pheno[, -1];

feats <- filt_feats;
pheno <- filt_pheno;

print(sprintf("Loaded %d pheno and %d feat data subjects, %d phenotypes.", nrow(feats), ncol(feats), length(pheno)));

feat_pheno <- cbind(feats, pheno);

print("Feature-phenotype column names:");
colnames(feat_pheno)

print("Null model parameters:");
model <- glm(pheno~., family="binomial", data=feat_pheno);
summary(model);

print("Writing phenotype predictions per null model.");
preds <- predict(model, feats, type="response");
write.table(preds, "R_preds.txt", sep="\t", quote=FALSE, col.names=FALSE);


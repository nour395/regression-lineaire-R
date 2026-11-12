# Régression linéaire : jeux de données A à J
# Lancer depuis la racine du dépôt (les données sont dans data/).
# Packages : MASS (fourni avec R), leaps, glmnet

set.seed(2026)
library(MASS)     # stepAIC
library(leaps)    # Cp de Mallows
library(glmnet)   # ridge et LASSO


# ----------------------------------------------------------------------
# Outils communs
# ----------------------------------------------------------------------

charger <- function(fichier) {
  nom <- load(file.path("data", fichier))
  get(nom)
}

# Vérifie les hypothèses sur le bruit à partir des résidus :
# espérance nulle, normalité (Q-Q plot + Kolmogorov-Smirnov sur les résidus studentisés),
# indépendance (résidus dans l'ordre) et variance constante (Scale-Location)
diagnostic_residus <- function(modele, titre = "") {
  n <- length(residuals(modele))
  ddl <- n - length(coef(modele)) - 1
  res_stud <- rstudent(modele)

  cat("Moyenne des résidus :", mean(residuals(modele)), "\n")
  print(ks.test(res_stud, "pt", df = ddl))

  op <- par(mfrow = c(2, 2))
  qqnorm(rstandard(modele), main = paste("Q-Q plot", titre))
  qqline(rstandard(modele), col = "red", lwd = 2)
  plot(ecdf(res_stud), do.points = FALSE, lwd = 2, main = "Kolmogorov : empirique vs Student")
  curve(pt(x, df = ddl), add = TRUE, col = "red", lwd = 2, lty = 2)
  plot(residuals(modele), main = "Résidus dans l'ordre des observations", ylab = "résidus")
  abline(h = 0, col = "red")
  plot(modele, which = 3)
  par(op)
}

# Nuage de points + modèle + intervalles de confiance et de prédiction à 95 %
trace_intervalles <- function(modele, df, titre = "") {
  plot(df$X, df$Y, pch = 16, col = "grey60", main = titre, xlab = "X", ylab = "Y")
  grille <- data.frame(X = seq(min(df$X), max(df$X), length.out = 300))
  ic <- predict(modele, newdata = grille, interval = "confidence")
  ip <- predict(modele, newdata = grille, interval = "prediction")
  lines(grille$X, ic[, "fit"], col = "red", lwd = 2)
  matlines(grille$X, ic[, c("lwr", "upr")], col = "blue", lty = 2, lwd = 2)
  matlines(grille$X, ip[, c("lwr", "upr")], col = "darkgreen", lty = 3, lwd = 2)
  legend("topleft", bty = "n", lty = c(1, 2, 3), lwd = 2, col = c("red", "blue", "darkgreen"),
         legend = c("modèle", "intervalle de confiance", "intervalle de prédiction"))
}


# ----------------------------------------------------------------------
# 1. Régression linéaire simple
# ----------------------------------------------------------------------

# --- Jeu A : le cas où tout va bien ---

A <- charger("data1.RData")
modele_A <- lm(Y ~ X, data = A)
summary(modele_A)

# Coefficients des moindres carrés recalculés à la main
x <- A$X; y <- A$Y; n <- length(x)
a_hat <- sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2)
b_hat <- mean(y) - a_hat * mean(x)
c(a_hat = a_hat, b_hat = b_hat)

# Résidus studentisés à la main : e_i / (sigma_(i) * sqrt(1 - h_ii))
X_mat <- cbind(1, x)
h <- diag(X_mat %*% solve(t(X_mat) %*% X_mat) %*% t(X_mat))
residus <- y - (b_hat + a_hat * x)
sigma_sans_i <- sapply(1:n, function(i) summary(lm(y[-i] ~ x[-i]))$sigma)
res_stud_main <- residus / (sigma_sans_i * sqrt(1 - h))
all.equal(unname(res_stud_main), unname(rstudent(modele_A)))

diagnostic_residus(modele_A, "jeu A")

# Intervalle de prédiction à la main (orange), superposé à celui de R
trace_intervalles(modele_A, A, "Jeu A")
x_new <- seq(min(x), max(x), length.out = 300)
sigma_hat <- sqrt(sum(residus^2) / (n - 2))
demi <- qt(0.975, n - 2) * sigma_hat * sqrt(1 + 1 / n + (x_new - mean(x))^2 / sum((x - mean(x))^2))
lines(x_new, b_hat + a_hat * x_new - demi, col = "orange", lty = 2, lwd = 3)
lines(x_new, b_hat + a_hat * x_new + demi, col = "orange", lty = 2, lwd = 3)

# --- Jeu B : bon R², mauvais bruit ---

B <- charger("data2.RData")
modele_B <- lm(Y ~ X, data = B)
summary(modele_B)$r.squared
diagnostic_residus(modele_B, "jeu B")

# --- Jeu C : la variance qui grandit ---

C <- charger("data3.RData")
modele_C <- lm(Y ~ X, data = C)
summary(modele_C)$r.squared
diagnostic_residus(modele_C, "jeu C")
trace_intervalles(modele_C, C, "Jeu C")

# --- Jeu D : une montée puis une descente, coupée en X = 32 ---

D <- charger("data4.RData")
plot(D$X, D$Y, main = "Jeu D", xlab = "X", ylab = "Y")

D_gauche <- D[D$X <= 32, ]
D_droite <- D[D$X > 32, ]
modele_g <- lm(Y ~ X, data = D_gauche)
modele_d <- lm(Y ~ X, data = D_droite)
c(R2_gauche = summary(modele_g)$r.squared, R2_droite = summary(modele_d)$r.squared)
diagnostic_residus(modele_g, "X <= 32")
diagnostic_residus(modele_d, "X > 32")

par(mfrow = c(1, 2))
trace_intervalles(modele_g, D_gauche, "X <= 32")
trace_intervalles(modele_d, D_droite, "X > 32")
par(mfrow = c(1, 1))


# ----------------------------------------------------------------------
# 2. Régression multiple et choix des variables
# ----------------------------------------------------------------------

# Compare Cp de Mallows, stepAIC, LASSO et ridge sur un jeu de données
selection_variables <- function(df, titre) {
  complet <- lm(Y ~ ., data = df)
  print(summary(complet))
  print(round(cor(df), 3))

  X <- as.matrix(df[, setdiff(names(df), "Y")])
  Y <- df$Y
  cp <- leaps(X, Y, method = "Cp")
  meilleurs <- order(cp$Cp)[1:5]
  cat("\nLes 5 meilleurs sous-modèles selon le Cp de Mallows :\n")
  print(cbind(cp$which[meilleurs, ], Cp = round(cp$Cp[meilleurs], 2)))

  vide <- lm(Y ~ 1, data = df)
  pas_a_pas <- stepAIC(vide, scope = list(lower = ~1, upper = formula(complet)), direction = "both", trace = FALSE)
  cat("\nModèle retenu par stepAIC :\n"); print(formula(pas_a_pas))

  lasso <- cv.glmnet(X, Y, alpha = 1)
  plot(lasso, main = paste("LASSO,", titre))
  cat("\nCoefficients du LASSO (lambda.1se) :\n"); print(coef(lasso, s = "lambda.1se"))

  # Ridge (alpha = 0) : les coefficients sont réduits mais aucun n'est mis à zéro
  ridge <- cv.glmnet(X, Y, alpha = 0)
  cat("\nCoefficients du ridge (lambda.1se) :\n"); print(coef(ridge, s = "lambda.1se"))
}

# --- Jeu E ---

E <- charger("data5.RData")
selection_variables(E, "jeu E")
summary(lm(Y ~ X1, data = E))$r.squared

# --- Jeu F ---

F_df <- charger("data6.RData")   # F est déjà le nom de FALSE en R
selection_variables(F_df, "jeu F")
summary(lm(Y ~ X1 + X3 + X4, data = F_df))$adj.r.squared

# --- Jeu G : des variables clones (X1 et X2) ---

G <- charger("data7.RData")
selection_variables(G, "jeu G")
summary(lm(Y ~ X5, data = G))$r.squared


# ----------------------------------------------------------------------
# 3. Variable qualitative : l'analyse de la variance
# ----------------------------------------------------------------------

anova_un_facteur <- function(df, titre) {
  df$X <- as.factor(df$X)
  boxplot(Y ~ X, data = df, col = "lightblue", main = titre, xlab = "groupe", ylab = "Y")
  points(seq_len(nlevels(df$X)), tapply(df$Y, df$X, mean), pch = 18, col = "red", cex = 1.6)
  modele <- lm(Y ~ X, data = df)
  print(anova(modele))                     # test de Fisher : toutes les moyennes égales ?
  print(summary(modele)$coefficients)      # écarts au premier groupe (référence)
  op <- par(mfrow = c(1, 2))
  qqnorm(residuals(modele)); qqline(residuals(modele), col = "red", lwd = 2)
  plot(fitted(modele), residuals(modele), xlab = "valeurs prédites", ylab = "résidus")
  abline(h = 0, col = "red")
  par(op)
}

# --- Jeu H ---

anova_un_facteur(charger("data8.RData"), "Jeu H")

# --- Jeu I ---

anova_un_facteur(charger("data9.RData"), "Jeu I")


# ----------------------------------------------------------------------
# 4. La régression polynomiale
# ----------------------------------------------------------------------

J <- charger("data10.RData")
plot(J$X, J$Y, pch = 16, col = "grey60", main = "Jeu J", xlab = "X", ylab = "Y")

# I() dit à R de calculer X^2 et X^3 tels quels
modele_J <- lm(Y ~ X + I(X^2) + I(X^3), data = J)
summary(modele_J)
diagnostic_residus(modele_J, "jeu J")
trace_intervalles(modele_J, J, "Jeu J : polynôme de degré 3")

# Test de Fisher : les termes X^2 et X^3 sont-ils utiles par rapport à une droite ?
anova(lm(Y ~ X, data = J), modele_J)

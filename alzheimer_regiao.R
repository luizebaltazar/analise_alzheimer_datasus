######################################################
# SÉRIE TEMPORAL - ALZHEIMER POR REGIÃO            #
# Brasil, 2010-2024                                 #
######################################################

# 1. PACOTES ----
library(forecast)
library(ggplot2)
library(tidyr)
library(dplyr)
library(zoo)
library(urca)
library(lmtest)
library(tseries)
library(prais)

# 2. IMPORTAR DADOS ----
dados <- read.csv("C:/Users/luize/Downloads/tabela6_alzheimer_regiao.csv")
dados$ano <- dados$Anos

# 3. GRÁFICO COM AS 5 REGIÕES JUNTAS ----
dados_long <- pivot_longer(dados,
                           cols = c(Norte, Nordeste, Sudeste, Sul, Centro_Oeste),
                           names_to = "regiao",
                           values_to = "taxa")

dados_long$regiao <- factor(dados_long$regiao,
                            levels = c("Norte", "Nordeste", "Sudeste",
                                       "Sul", "Centro_Oeste"),
                            labels = c("Norte", "Nordeste", "Sudeste",
                                       "Sul", "Centro-Oeste"))

ggplot(dados_long, aes(x = Anos, y = taxa, color = regiao)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Norte"        = "#E41A1C",
                                "Nordeste"     = "#FF7F00",
                                "Sudeste"      = "#4DAF4A",
                                "Sul"          = "#377EB8",
                                "Centro-Oeste" = "#984EA3")) +
  facet_wrap(~ regiao, scales = "free_y", ncol = 1) +
  labs(title = "Taxa de Mortalidade por Doença de Alzheimer",
       subtitle = "Segundo região — Brasil, 2010–2024",
       x = "Ano",
       y = "Taxa por 100.000 habitantes",
       caption = "Fonte: SIM/DATASUS e IBGE, Projeção Populacional, 2024") +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

# ============================================================
# 4. ANÁLISE - NORTE ----
# ============================================================
dados$tx  <- dados$Norte
dados$log <- log(dados$tx)

st_n    <- ts(data = dados$tx,  start = c(2010), frequency = 1)
logst_n <- ts(data = dados$log, start = c(2010), frequency = 1)

plot(st_n,
     main = "Taxa de Mortalidade por Alzheimer - Norte",
     xlab = "Ano", ylab = "Taxa por 100.000 habitantes", type = "o")

m_n <- prais_winsten(log ~ ano, index = "ano", data = dados)
summary(m_n)
acf(diff(logst_n),  main = "FAC - Norte")
Box.test(residuals(m_n), lag = 1, type = "Ljung-Box")

# ============================================================
# 5. ANÁLISE - NORDESTE ----
# ============================================================
dados$tx  <- dados$Nordeste
dados$log <- log(dados$tx)

st_ne    <- ts(data = dados$tx,  start = c(2010), frequency = 1)
logst_ne <- ts(data = dados$log, start = c(2010), frequency = 1)

plot(st_ne,
     main = "Taxa de Mortalidade por Alzheimer - Nordeste",
     xlab = "Ano", ylab = "Taxa por 100.000 habitantes", type = "o")

m_ne <- prais_winsten(log ~ ano, index = "ano", data = dados)
summary(m_ne)
acf(diff(logst_ne), main = "FAC - Nordeste")
Box.test(residuals(m_ne), lag = 1, type = "Ljung-Box")

# ============================================================
# 6. ANÁLISE - SUDESTE ----
# ============================================================
dados$tx  <- dados$Sudeste
dados$log <- log(dados$tx)

st_se    <- ts(data = dados$tx,  start = c(2010), frequency = 1)
logst_se <- ts(data = dados$log, start = c(2010), frequency = 1)

plot(st_se,
     main = "Taxa de Mortalidade por Alzheimer - Sudeste",
     xlab = "Ano", ylab = "Taxa por 100.000 habitantes", type = "o")

m_se <- prais_winsten(log ~ ano, index = "ano", data = dados)
summary(m_se)
acf(diff(logst_se), main = "FAC - Sudeste")
Box.test(residuals(m_se), lag = 1, type = "Ljung-Box")

# ============================================================
# 7. ANÁLISE - SUL ----
# ============================================================
dados$tx  <- dados$Sul
dados$log <- log(dados$tx)

st_s    <- ts(data = dados$tx,  start = c(2010), frequency = 1)
logst_s <- ts(data = dados$log, start = c(2010), frequency = 1)

plot(st_s,
     main = "Taxa de Mortalidade por Alzheimer - Sul",
     xlab = "Ano", ylab = "Taxa por 100.000 habitantes", type = "o")

m_s <- prais_winsten(log ~ ano, index = "ano", data = dados)
summary(m_s)
acf(diff(logst_s),  main = "FAC - Sul")
Box.test(residuals(m_s), lag = 1, type = "Ljung-Box")

# ============================================================
# 8. ANÁLISE - CENTRO-OESTE ----
# ============================================================
dados$tx  <- dados$Centro_Oeste
dados$log <- log(dados$tx)

st_co    <- ts(data = dados$tx,  start = c(2010), frequency = 1)
logst_co <- ts(data = dados$log, start = c(2010), frequency = 1)

plot(st_co,
     main = "Taxa de Mortalidade por Alzheimer - Centro-Oeste",
     xlab = "Ano", ylab = "Taxa por 100.000 habitantes", type = "o")

m_co <- prais_winsten(log ~ ano, index = "ano", data = dados)
summary(m_co)
acf(diff(logst_co), main = "FAC - Centro-Oeste")
Box.test(residuals(m_co), lag = 1, type = "Ljung-Box")

# ============================================================
# 9. CALCULAR VPA E IC95% ----
# ============================================================
calcular_vpa <- function(modelo, regiao) {
  b1     <- coef(modelo)["ano"]
  se_b1  <- summary(modelo)$coefficients["ano", "Std. Error"]
  t_crit <- qt(0.975, df = 13)

  ic_inf <- b1 - t_crit * se_b1
  ic_sup <- b1 + t_crit * se_b1

  vpa     <- (exp(b1)     - 1) * 100
  vpa_inf <- (exp(ic_inf) - 1) * 100
  vpa_sup <- (exp(ic_sup) - 1) * 100
  p_valor <- summary(modelo)$coefficients["ano", "Pr(>|t|)"]

  cat("\n", rep("=", 50), "\n", sep = "")
  cat(" Região:", regiao, "\n")
  cat(rep("=", 50), "\n", sep = "")
  cat(sprintf(" VPA    = %+.2f%% ao ano\n", vpa))
  cat(sprintf(" IC95%%  = %+.2f%% a %+.2f%%\n", vpa_inf, vpa_sup))
  cat(sprintf(" p-valor = %.4f\n", p_valor))

  return(data.frame(Regiao  = regiao,
                    VPA     = round(vpa, 2),
                    IC_inf  = round(vpa_inf, 2),
                    IC_sup  = round(vpa_sup, 2),
                    p_valor = round(p_valor, 4)))
}

res_n  <- calcular_vpa(m_n,  "Norte")
res_ne <- calcular_vpa(m_ne, "Nordeste")
res_se <- calcular_vpa(m_se, "Sudeste")
res_s  <- calcular_vpa(m_s,  "Sul")
res_co <- calcular_vpa(m_co, "Centro-Oeste")

# 10. TABELA RESUMO ----
tabela_vpa_regiao <- rbind(res_n, res_ne, res_se, res_s, res_co)
cat("\n\n========== TABELA RESUMO ==========\n")
print(tabela_vpa_regiao, row.names = FALSE)

# Exportar tabela
write.csv(tabela_vpa_regiao, "C:/Users/luize/Downloads/tabela_vpa_alzheimer_regiao.csv",
          row.names = FALSE)
cat("\n✓ Tabela exportada: tabela_vpa_alzheimer_regiao.csv\n")

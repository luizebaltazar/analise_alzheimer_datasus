######################################################
# SÉRIE TEMPORAL - ALZHEIMER POR FAIXA ETÁRIA      #
# Brasil, 2010-2024                                 #
######################################################

# 1. PACOTES ----
library(forecast)
library(ggplot2)
library(tidyr)
library(zoo)
library(urca)
library(lmtest)
library(tseries)
library(prais)

# 2. IMPORTAR DADOS ----
dados <- read.csv("C:/Users/luize/Downloads/tabela2_alzheimer_faixa_etaria.csv")

# 3. GRÁFICO COM AS 3 FAIXAS JUNTAS ----
dados_long <- pivot_longer(dados,
                           cols = c(fx60_69, fx70_79, fx80mais),
                           names_to = "faixa",
                           values_to = "taxa")

dados_long$faixa <- factor(dados_long$faixa,
                           levels = c("fx60_69", "fx70_79", "fx80mais"),
                           labels = c("60-69 anos", "70-79 anos", "≥80 anos"))

ggplot(dados_long, aes(x = Anos, y = taxa, color = faixa)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("steelblue", "tomato", "darkgreen")) +
  facet_wrap(~ faixa, scales = "free_y", ncol = 1) +
  labs(title = "Taxa de Mortalidade por Doença de Alzheimer",
       subtitle = "Segundo faixa etária — Brasil, 2010–2024",
       x = "Ano",
       y = "Taxa por 100.000 habitantes",
       caption = "Fonte: Sistema de Informações sobre Mortalidade (SIM/DATASUS) e IBGE, Projeção Populacional, 2024") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

# ============================================================
# 4. ANÁLISE POR FAIXA ETÁRIA - 60 A 69 ANOS ----
# ============================================================
dados$tx  <- dados$fx60_69
dados$log <- log(dados$tx)
dados$ano <- dados$Anos

st_60    <- ts(data = dados$tx,  start = c(2010), frequency = 1)
logst_60 <- ts(data = dados$log, start = c(2010), frequency = 1)

plot(st_60,
     main = "Taxa de Mortalidade por Alzheimer - 60 a 69 anos",
     xlab = "Ano",
     ylab = "Taxa por 100.000 habitantes",
     type = "o")

autoplot(logst_60) +
  labs(title = "Série em Log - 60 a 69 anos", x = "Ano", y = "Log da Taxa")

m_60 <- prais_winsten(log ~ ano, index = "ano", data = dados)
summary(m_60)

acf(diff(logst_60), main = "FAC - 60 a 69 anos")
acf(diff(logst_60), main = "FACP - 60 a 69 anos", type = "partial")
Box.test(residuals(m_60), lag = 1, type = "Ljung-Box")

# ============================================================
# 5. ANÁLISE POR FAIXA ETÁRIA - 70 A 79 ANOS ----
# ============================================================
dados$tx  <- dados$fx70_79
dados$log <- log(dados$tx)

st_70    <- ts(data = dados$tx,  start = c(2010), frequency = 1)
logst_70 <- ts(data = dados$log, start = c(2010), frequency = 1)

plot(st_70,
     main = "Taxa de Mortalidade por Alzheimer - 70 a 79 anos",
     xlab = "Ano",
     ylab = "Taxa por 100.000 habitantes",
     type = "o")

autoplot(logst_70) +
  labs(title = "Série em Log - 70 a 79 anos", x = "Ano", y = "Log da Taxa")

m_70 <- prais_winsten(log ~ ano, index = "ano", data = dados)
summary(m_70)

acf(diff(logst_70), main = "FAC - 70 a 79 anos")
acf(diff(logst_70), main = "FACP - 70 a 79 anos", type = "partial")
Box.test(residuals(m_70), lag = 1, type = "Ljung-Box")

# ============================================================
# 6. ANÁLISE POR FAIXA ETÁRIA - 80 ANOS OU MAIS ----
# ============================================================
dados$tx  <- dados$fx80mais
dados$log <- log(dados$tx)

st_80    <- ts(data = dados$tx,  start = c(2010), frequency = 1)
logst_80 <- ts(data = dados$log, start = c(2010), frequency = 1)

plot(st_80,
     main = "Taxa de Mortalidade por Alzheimer - ≥80 anos",
     xlab = "Ano",
     ylab = "Taxa por 100.000 habitantes",
     type = "o")

autoplot(logst_80) +
  labs(title = "Série em Log - ≥80 anos", x = "Ano", y = "Log da Taxa")

m_80 <- prais_winsten(log ~ ano, index = "ano", data = dados)
summary(m_80)

acf(diff(logst_80), main = "FAC - ≥80 anos")
acf(diff(logst_80), main = "FACP - ≥80 anos", type = "partial")
Box.test(residuals(m_80), lag = 1, type = "Ljung-Box")

# ============================================================
# 7. CALCULAR VPA E IC95% PARA CADA FAIXA ----
# ============================================================
calcular_vpa <- function(modelo, faixa) {
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
  cat(" Faixa:", faixa, "\n")
  cat(rep("=", 50), "\n", sep = "")
  cat(sprintf(" VPA    = %+.2f%% ao ano\n", vpa))
  cat(sprintf(" IC95%%  = %+.2f%% a %+.2f%%\n", vpa_inf, vpa_sup))
  cat(sprintf(" p-valor = %.4f\n", p_valor))
  
  return(data.frame(Faixa = faixa,
                    VPA = round(vpa, 2),
                    IC_inf = round(vpa_inf, 2),
                    IC_sup = round(vpa_sup, 2),
                    p_valor = round(p_valor, 4)))
}

res_60 <- calcular_vpa(m_60, "60-69 anos")
res_70 <- calcular_vpa(m_70, "70-79 anos")
res_80 <- calcular_vpa(m_80, "≥80 anos")

# 8. TABELA RESUMO ----
tabela_vpa <- rbind(res_60, res_70, res_80)
cat("\n\n========== TABELA RESUMO ==========\n")
print(tabela_vpa, row.names = FALSE)

# Exportar tabela
write.csv(tabela_vpa, "C:/Users/luize/Downloads/tabela_vpa_alzheimer.csv",
          row.names = FALSE)
cat("\n✓ Tabela exportada: tabela_vpa_alzheimer.csv\n")


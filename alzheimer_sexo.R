######################################################
# ROTEIRO BÁSICO PARA REALIZAR SÉRIES TEMPORAIS ######
######################################################
# 1.Pacotes necessários para as séries temporais
#pacotes essencias utilizados
install.packages("forecast")
install.packages("ggplot2")
install.packages("urca")
install.packages("zoo")
install.packages("lmtest")
install.packages("seasonal")
install.packages("tseries")
install.packages("prais")
install.packages("timeSeries")
### depois de instalados, ativar os pacotes ###
library(forecast)
library(ggplot2)
library(seasonal)
library(timeSeries)
library(zoo)
library(urca)
library(lmtest)
library(tseries)
library(prais)
######################################################
# SÉRIE TEMPORAL - ALZHEIMER POR SEXO              #
# Brasil, 2010-2024                                 #
######################################################


# 2. IMPORTAR DADOS ----
dados <- read.csv("C:/Users/luize/Downloads/tabela4_alzheimer_sexo.csv")
dados$ano <- dados$Anos

# 3. GRÁFICO COM OS 2 SEXOS JUNTOS ----
dados_long <- pivot_longer(dados,
                           cols = c(masculino, feminino),
                           names_to = "sexo",
                           values_to = "taxa")

dados_long$sexo <- factor(dados_long$sexo,
                          levels = c("masculino", "feminino"),
                          labels = c("Masculino", "Feminino"))

ggplot(dados_long, aes(x = Anos, y = taxa, color = sexo)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#2166AC", "#D6604D")) +
  facet_wrap(~ sexo, scales = "free_y", ncol = 1) +
  labs(title = "Taxa de Mortalidade por Doença de Alzheimer",
       subtitle = "Segundo sexo — Brasil, 2010–2024",
       x = "Ano",
       y = "Taxa por 100.000 habitantes",
       caption = "Fonte: Sistema de Informações sobre Mortalidade (SIM/DATASUS) e IBGE, Projeção Populacional, 2024") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

# ============================================================
# 4. ANÁLISE - MASCULINO ----
# ============================================================
dados$tx  <- dados$masculino
dados$log <- log(dados$tx)

st_m    <- ts(data = dados$tx,  start = c(2010), frequency = 1)
logst_m <- ts(data = dados$log, start = c(2010), frequency = 1)

plot(st_m,
     main = "Taxa de Mortalidade por Alzheimer - Masculino",
     xlab = "Ano",
     ylab = "Taxa por 100.000 habitantes",
     type = "o")

autoplot(logst_m) +
  labs(title = "Série em Log - Masculino", x = "Ano", y = "Log da Taxa")

m_masc <- prais_winsten(log ~ ano, index = "ano", data = dados)
summary(m_masc)

acf(diff(logst_m),  main = "FAC - Masculino")
acf(diff(logst_m),  main = "FACP - Masculino", type = "partial")
Box.test(residuals(m_masc), lag = 1, type = "Ljung-Box")

# ============================================================
# 5. ANÁLISE - FEMININO ----
# ============================================================
dados$tx  <- dados$feminino
dados$log <- log(dados$tx)

st_f    <- ts(data = dados$tx,  start = c(2010), frequency = 1)
logst_f <- ts(data = dados$log, start = c(2010), frequency = 1)

plot(st_f,
     main = "Taxa de Mortalidade por Alzheimer - Feminino",
     xlab = "Ano",
     ylab = "Taxa por 100.000 habitantes",
     type = "o")

autoplot(logst_f) +
  labs(title = "Série em Log - Feminino", x = "Ano", y = "Log da Taxa")

m_fem <- prais_winsten(log ~ ano, index = "ano", data = dados)
summary(m_fem)

acf(diff(logst_f),  main = "FAC - Feminino")
acf(diff(logst_f),  main = "FACP - Feminino", type = "partial")
Box.test(residuals(m_fem), lag = 1, type = "Ljung-Box")

# ============================================================
# RAZÃO FEMININO/MASCULINO POR ANO ----
# ============================================================

dados$razao <- dados$feminino / dados$masculino

# Imprimir
cat("\n========== RAZÃO FEMININO/MASCULINO ==========\n")
print(data.frame(Ano = dados$Anos, Razao = round(dados$razao, 2)), row.names = FALSE)

# Gráfico
ggplot(dados, aes(x = Anos, y = razao)) +
  geom_line(linewidth = 1, color = "purple") +
  geom_point(size = 2, color = "purple") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  labs(title = "Razão de Mortalidade por Alzheimer: Feminino/Masculino",
       subtitle = "Brasil, 2010–2024",
       x = "Ano",
       y = "Razão (Feminino/Masculino)",
       caption = "Fonte: SIM/DATASUS e IBGE, Projeção Populacional, 2024") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

# ============================================================
# 6. CALCULAR VPA E IC95% ----
# ============================================================
calcular_vpa <- function(modelo, grupo) {
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
  cat(" Sexo:", grupo, "\n")
  cat(rep("=", 50), "\n", sep = "")
  cat(sprintf(" VPA    = %+.2f%% ao ano\n", vpa))
  cat(sprintf(" IC95%%  = %+.2f%% a %+.2f%%\n", vpa_inf, vpa_sup))
  cat(sprintf(" p-valor = %.4f\n", p_valor))
  
  return(data.frame(Sexo = grupo,
                    VPA = round(vpa, 2),
                    IC_inf = round(vpa_inf, 2),
                    IC_sup = round(vpa_sup, 2),
                    p_valor = round(p_valor, 4)))
}

res_masc <- calcular_vpa(m_masc, "Masculino")
res_fem  <- calcular_vpa(m_fem,  "Feminino")

# 7. TABELA RESUMO ----
tabela_vpa_sexo <- rbind(res_masc, res_fem)
cat("\n\n========== TABELA RESUMO ==========\n")
print(tabela_vpa_sexo, row.names = FALSE)

# Exportar tabela
write.csv(tabela_vpa_sexo, "C:/Users/luize/Downloads/tabela_vpa_alzheimer_sexo.csv",
          row.names = FALSE)
cat("\n✓ Tabela exportada: tabela_vpa_alzheimer_sexo.csv\n")

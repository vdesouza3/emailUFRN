# Carregar pacotes necessários
library(ncdf4)
library(ggplot2)
library(dplyr)
library(terra)
library(sf)

# --- Configurações Iniciais ---
# Caminho para os arquivos e diretórios
caminho <- "/home/laico/Área de Trabalho/VRSS/dados selecionados"
dir_graficos <- "/home/laico/Área de Trabalho/VRSS/graficos_linha_selecionados/"
dir.create(dir_graficos, showWarnings = FALSE, recursive = TRUE)

# Definir anos, alturas e os estados de interesse
anos_desejados <- c(2003)
alturas <- c(120, 150, 180)
# CORREÇÃO: Usando a grafia com acento para Goias
estados_desejados <- c("RIO GRANDE DO SUL") 

# Carregar shapefile e converter os nomes dos estados para maiúsculas
shapefile_estados <- st_read("/home/laico/Área de Trabalho/VRSS/documentos/BR_UF_2024/BR_UF_2024.shp", quiet = TRUE) |>
  st_transform(crs = 4326) |>
  mutate(NM_UF = toupper(NM_UF))

# --- Funções de Processamento ---

# Função para interpolar u e v
interpolar_uv <- function(u10, v10, u100, v100, z) {
  u_z <- u10 + (log(z/10) / log(100/10)) * (u100 - u10)
  v_z <- v10 + (log(z/10) / log(100/10)) * (v100 - v10)
  list(u = u_z, v = v_z)
}

# Função para processar UM MÊS e retornar os dados
processar_mes_selecionado <- function(mes) {
  message(sprintf("🔄 Processando dados para o mês de %s...", mes))
  
  arquivo_mes <- list.files(caminho, pattern = sprintf("^%s.*\\.nc$", mes), full.names = TRUE)
  
  if (length(arquivo_mes) == 0) {
    message(sprintf("⚠️ Nenhum arquivo encontrado para o mês de %s. Pulando...", mes))
    return(NULL)
  }
  
  nc <- nc_open(arquivo_mes)
  anos_disponiveis <- 1961:2024
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  todos_dados <- list()
  
  for (i in which(anos_disponiveis %in% anos_desejados)) {
    ano_atual <- anos_disponiveis[i]
    
    u10 <- ncvar_get(nc, "10u", start = c(1, 1, i), count = c(-1, -1, 1))
    v10 <- ncvar_get(nc, "10v", start = c(1, 1, i), count = c(-1, -1, 1))
    u100 <- ncvar_get(nc, "100u", start = c(1, 1, i), count = c(-1, -1, 1))
    v100 <- ncvar_get(nc, "100v", start = c(1, 1, i), count = c(-1, -1, 1))
    
    for (h in alturas) {
      interp <- interpolar_uv(u10, v10, u100, v100, h)
      vel <- sqrt(interp$u^2 + interp$v^2)
      
      df_temp <- expand.grid(lon = lon, lat = lat) |>
        mutate(
          vel = as.vector(vel),
          ano = ano_atual,
          mes = mes,  
          altura = h
        ) |>
        filter(!is.na(vel))
      
      df_sf <- st_as_sf(df_temp, coords = c("lon", "lat"), crs = 4326)
      df_join <- st_join(df_sf, shapefile_estados, join = st_intersects) |>
        as.data.frame() |>
        filter(NM_UF %in% toupper(estados_desejados))
      
      if(nrow(df_join) > 0) {
        todos_dados[[length(todos_dados) + 1]] <- df_join
      }
    }
  }
  nc_close(nc)
  return(bind_rows(todos_dados))
}

# --- Execução do Script ---
message("Script iniciado em: ", Sys.time())
message("⚠️ Os dados do arquivo parecem ser médias mensais. O gráfico mostrará a média do mês, não a variação diária.")
message("🔄 Processando todos os meses do ano de 2000...")

meses_do_ano <- c("jan", "fev", "mar", "abr", "mai", "jun", "jul", "ago", "set", "out", "nov", "dez")
dados_de_todos_os_meses <- list()

for (m in meses_do_ano) {
  dados_mes <- processar_mes_selecionado(m)
  if (!is.null(dados_mes)) {
    dados_de_todos_os_meses[[m]] <- dados_mes
  }
}

# Combina todos os dados de todos os meses em um único dataframe
dados_finais <- bind_rows(dados_de_todos_os_meses)

if (nrow(dados_finais) == 0) {
  message("⚠️ Nenhum dado encontrado para os estados e meses selecionados.")
} else {
  # Agrega os dados pela altura, mês e estado
  dados_agregados <- dados_finais |>
    group_by(mes, altura, NM_UF) |>
    summarise(vel_media = mean(vel, na.rm = TRUE), .groups = "drop") |>
    mutate(altura = as.factor(altura),
           mes = factor(mes, levels = meses_do_ano)) 
  
  # Cria o gráfico de linha com painéis para cada estado
  p <- ggplot(dados_agregados, aes(x = mes, y = vel_media, color = altura, group = altura)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    facet_wrap(~NM_UF, scales = "free_y") +
    labs(
      title = sprintf("Velocidade Média Mensal do Vento - %d", anos_desejados[1]),
      subtitle = "Comparação entre diferentes alturas e meses",
      x = "Mês",
      y = "Velocidade do Vento (m/s)",
      color = "Altura (m)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )
  
  # Converte os nomes dos estados para uma única string para o nome do arquivo
  regioes_string <- paste(estados_desejados, collapse = "-")
  
  nome_img <- sprintf("vento_media_mensal_%s_%d.png", regioes_string, anos_desejados[1])
  
  ggsave(filename = file.path(dir_graficos, nome_img), plot = p, width = 12, height = 9)
  
  cat(sprintf("✅ Gráfico salvo: %s\n", nome_img))
}

message("\n✅ Processamento completo em: ", Sys.time())



############################################################################################


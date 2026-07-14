#Script CC e CCET#

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from google.colab import files

# ====================================================================
# PARTE 1: CARREGAMENTO E LEITURA
# ====================================================================

print(">>> Por favor, carregue o arquivo do Período SECO:")
uploaded_seco = files.upload()
nome_arquivo_seco = list(uploaded_seco.keys())[0]

print("\n>>> Por favor, carregue o arquivo do Período CHUVOSO:")
uploaded_chuvoso = files.upload()
nome_arquivo_chuvoso = list(uploaded_chuvoso.keys())[0]

try:
    df_seco = pd.read_excel(nome_arquivo_seco)
    df_chuvoso = pd.read_excel(nome_arquivo_chuvoso)
except Exception:
    df_seco = pd.read_csv(nome_arquivo_seco, encoding='latin-1', sep=None, engine='python')
    df_chuvoso = pd.read_csv(nome_arquivo_chuvoso, encoding='latin-1', sep=None, engine='python')

df_seco['Período'] = 'Seco'
df_chuvoso['Período'] = 'Chuvoso'
df_combinado = pd.concat([df_seco, df_chuvoso], ignore_index=True)
df_combinado.columns = df_combinado.columns.str.strip()

# ====================================================================
# PARTE 2: RECÁLCULO DO TEV E AGRUPAMENTO (MÉDIAS)
# ====================================================================

colunas_numericas = ['Temperatura do ar (oC)', 'Umidade Relativa (%)', 'Velocidade do Vento (m/s)', 
                     'ITU (oC)', 'IDT (oC)', 'TEV (oC)']

for col in colunas_numericas:
    df_combinado[col] = pd.to_numeric(df_combinado[col], errors='coerce')

# Recálculo do TEV para precisão
Ta = df_combinado['Temperatura do ar (oC)']
UR = df_combinado['Umidade Relativa (%)']
Vvel = df_combinado['Velocidade do Vento (m/s)']
denominador = 0.68 - (0.0014 * UR) + (1 / (1.76 + 1.4 * (Vvel**0.75)))
df_combinado['TEV (oC)'] = 37 - ((37 - Ta) / denominador) - (0.29 * Ta * (1 - (UR / 100)))

# Mapeamento e Agrupamento
mapeamento_local = {'TH1': 'CCUF', 'TH2': 'CCUF', 'TH3': 'CCET', 'TH4': 'CCET'}
df_combinado['Grupo_Medida'] = df_combinado['Termohigrômetro'].map(mapeamento_local)

df_medias = df_combinado.groupby(
    ['Dias', 'Mês', 'Ano', 'Hora_min', 'Período', 'Grupo_Medida']
)[colunas_numericas].mean().reset_index()

# ====================================================================
# PARTE 3: CONFIGURAÇÕES E FUNÇÃO DE PLOTAGEM
# ====================================================================

limites_conforto = {
    'ITU (oC)': {'min': 21.0, 'max': 24.0},
    'IDT (oC)': {'max': 24.0},
    'TEV (oC)': {'min': 22.0, 'max': 25.0}
}

meanprops = {'marker': 'o', 'markeredgecolor': 'black', 'markerfacecolor': 'red', 'markersize': 8}
paleta = {'Seco': 'darkorange', 'Chuvoso': 'skyblue'}

def criar_painel(lista_variaveis, nome_saida):
    fig, axes = plt.subplots(1, 3, figsize=(24, 6))
    letras = ['(A)', '(B)', '(C)']
    
    proxy_mean = plt.plot([], [], marker='o', markeredgecolor='black', 
                          markerfacecolor='red', markersize=8, linestyle='None', label='Média')[0]

    for i, variavel in enumerate(lista_variaveis):
        ax = axes[i]
        sns.boxplot(ax=ax, x='Grupo_Medida', y=variavel, hue='Período', data=df_medias,
                    order=['CCUF', 'CCET'], palette=paleta, width=0.4, gap=0.1,
                    medianprops={'color': 'black', 'linewidth': 2},
                    showmeans=True, meanprops=meanprops)
        
        # Título lateral/Letra
        ax.set_title(letras[i], loc='left', fontweight='bold', fontsize=14)
        ax.set_xlabel('Local Agrupado (CCUF e CCET)', fontsize=11)
        ax.set_ylabel(variavel, fontsize=11)
        ax.grid(axis='y', linestyle='--', alpha=0.6)

        # Adição de Limites de Conforto
        handles_conforto = []
        labels_conforto = []
        if variavel in limites_conforto:
            lim = limites_conforto[variavel]
            if 'min' in lim:
                l1 = ax.axhline(lim['min'], color='red', ls='--', lw=1.5)
                handles_conforto.append(l1); labels_conforto.append('Limite Superior')
            if 'max' in lim:
                label_max = 'Limite (< 24°C)' if variavel == 'IDT (oC)' else 'Limite Inferior'
                l2 = ax.axhline(lim['max'], color='red', ls='--', lw=1.5)
                handles_conforto.append(l2); labels_conforto.append(label_max)

# --- Gerenciamento da Legenda Apenas na ÚLTIMA imagem (i == 2) ---
        if i == 2:
            # Pegamos os handles e labels que o gráfico gerou automaticamente
            h_raw, l_raw = ax.get_legend_handles_labels()
            
            # Criamos dicionário para remover duplicatas mantendo a ordem
            # Isso remove o 'Média' extra gerado pelo showmeans=True
            dict_legenda = dict(zip(l_raw, h_raw))
            
            # Montamos a lista final manualmente para garantir a ordem correta
            all_h = [proxy_mean]
            all_l = ['Média']
            
            # Adiciona Períodos (Seco e Chuvoso)
            for p in ['Seco', 'Chuvoso']:
                if p in dict_legenda:
                    all_h.append(dict_legenda[p])
                    all_l.append(p)
            
            # Adiciona os Limites que foram capturados durante o loop
            all_h.extend(handles_conforto)
            all_l.extend(labels_conforto)
            
            ax.legend(all_h, all_l, title='Legenda', 
                      bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
        else:
            if ax.get_legend():
                ax.get_legend().remove()

    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(nome_saida, dpi=300)
    plt.show()
    files.download(nome_saida)

# ====================================================================
# PARTE 4: EXECUÇÃO DOS PAINÉIS
# ====================================================================

# Painel 1: Variáveis Meteorológicas
vars_clima = ['Temperatura do ar (oC)', 'Umidade Relativa (%)', 'Velocidade do Vento (m/s)']
print("\n>>> Gerando Painel de Variáveis Climáticas...")
criar_painel(vars_clima, 'Painel_Variaveis_Climaticas.png')

# Painel 2: Índices de Conforto
vars_indices = ['ITU (oC)', 'IDT (oC)', 'TEV (oC)']
print("\n>>> Gerando Painel de Índices de Conforto...")
criar_painel(vars_indices, 'Painel_Indices_Conforto.png')

print("\nConcluído! Os dois painéis foram gerados e o download iniciado.")

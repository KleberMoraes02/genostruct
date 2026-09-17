import streamlit as st
import requests
import pandas as pd
import py3Dmol
from stmol import showmol
import os

# --- CONFIGURAÇÃO GLOBAL DA PÁGINA ---
st.set_page_config(
    page_title="GenoStruct MVP",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- DICIONÁRIO DE APOIO ---
MAPA_AMINOACIDOS = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
}

# --- FUNÇÃO DE DADOS ---
@st.cache_data(show_spinner=False)
def carregar_banco_mutacoes():
    try:
        # Lê o arquivo direto da pasta 'dados' que criamos
        caminho = os.path.join("dados", "banco_teste.csv.gz")
        return pd.read_csv(caminho)
    except Exception:
        return pd.DataFrame()

@st.cache_resource(show_spinner=False)
def buscar_pdb_alphafold(uniprot_id):
    try:
        url_api = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
        resposta_api = requests.get(url_api, timeout=5)
        if resposta_api.status_code == 200 and len(resposta_api.json()) > 0:
            url_pdb = resposta_api.json()[0]['pdbUrl']
            resposta_pdb = requests.get(url_pdb, timeout=5)
            if resposta_pdb.status_code == 200:
                return resposta_pdb.text
        return None
    except:
        return None

# --- MENU LATERAL (SIDEBAR) ---
with st.sidebar:
    st.image("https://cdn-icons-png.flaticon.com/512/3022/3022421.png", width=60)
    st.title("GenoStruct")
    st.caption("v1.1.0 (Filtro de Colunas)")
    st.markdown("---")
    st.markdown("**⚙️ Credenciais de Modelagem**")
    token_swiss = st.text_input("SWISS-MODEL API Token:", type="password")

# --- CORPO PRINCIPAL DO SITE ---
st.title("GenoStruct: Integração Genômica e Estrutural")
st.markdown("Busque um gene para correlacionar variantes clínicas com predições estruturais.")

df_mutacoes = carregar_banco_mutacoes()

if not df_mutacoes.empty:
    gene_buscado_raw = st.text_input("Nome do Gene alvo (ex: ABCD2):")
    gene_buscado = gene_buscado_raw.strip().upper()

    if gene_buscado:
        if len(gene_buscado) < 2:
            st.warning("⚠️ Digite um nome de gene válido com pelo menos 2 letras.")
        else:
            with st.spinner("Sincronizando com UniProtKB..."):
                try:
                    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_buscado} AND organism_id:9606 AND reviewed:true&format=json&size=1"
                    resposta = requests.get(url, timeout=10) 

                    if resposta.status_code == 200 and len(resposta.json()['results']) > 0:
                        proteina = resposta.json()['results'][0]
                        uniprot_id = proteina['primaryAccession']
                        sequencia_selvagem = proteina['sequence']['value']

                        st.success(f"Conexão estabelecida! Alvo: **{gene_buscado}** (Accession: {uniprot_id})")

                        # --- 1. PAINEL CLÍNICO (TELA CHEIA) ---
                        st.subheader("📊 Perfil Mutacional Clínico")
                        mutacoes_filtradas = df_mutacoes[df_mutacoes['Gene'].str.upper() == gene_buscado]

                        # --- NOVO: DICIONÁRIO DE COLUNAS (Sanfona) ---
                        with st.expander("📖 Dicionário de Colunas (Clique para expandir)"):
                            st.markdown("""
                            Este guia explica os dados genômicos, clínicos e preditivos exibidos na tabela abaixo.
                            
                            **Identificação e Genética:**
                            - **Gene:** Símbolo oficial aprovado pelo comitê HGNC.
                            - **HGNC_ID:** Identificador numérico único do gene.
                            - **UniProt_Accession:** Código do gene no banco mundial de proteínas (UniProt).
                            - **rsID:** Identificador universal da variante no banco dbSNP (Reference SNP cluster ID).
                            
                            **Nomenclatura (Padrão HGVS):**
                            - **HGVS_c:** Alteração no nível do DNA (sequência codificante). *Ex: c.2219C>T (Citosina trocada por Timina).*
                            - **HGVS_p / Variante:** Alteração no nível da Proteína. *Ex: p.Ser740Phe (Serina trocada por Fenilalanina na posição 740).*
                            - **Posicao:** A posição numérica exata do aminoácido afetado na cadeia proteica.
                            
                            **Predição e Relevância Clínica:**
                            - **AlphaMissense_Class:** Inteligência artificial do Google DeepMind que prevê se a mutação é Patogênica, Benigna ou Ambígua.
                            - **ClinVar:** Registro de relevância clínica real observada em pacientes médicos.
                            - **Frequency (gnomAD):** A frequência com que esta mutação ocorre na população mundial geral.
                            """)
                        # ---------------------------------------------

                        if not mutacoes_filtradas.empty:
                            
                            # --- NOVO: SELETOR DE COLUNAS ---
                            todas_as_colunas = mutacoes_filtradas.columns.tolist()
                            
                            # Escolhemos algumas colunas para virem selecionadas por padrão
                            colunas_padrao = ['Gene', 'Variante', 'Posicao', 'Troca', 'AlphaMissense_Class', 'ClinVar']
                            colunas_padrao = [c for c in colunas_padrao if c in todas_as_colunas] # Evita erros se a coluna não existir

                            colunas_selecionadas = st.multiselect(
                                "Selecione as colunas que deseja visualizar na tabela:",
                                options=todas_as_colunas,
                                default=colunas_padrao
                            )

                            # A tabela agora só mostra o que o usuário selecionou no menu acima
                            df_exibicao = mutacoes_filtradas[colunas_selecionadas]
                            st.dataframe(df_exibicao, use_container_width=True, hide_index=True)
                            # --------------------------------

                            st.markdown("---")
                            
                            # --- 2. PAINEL DE SIMULAÇÃO E 3D (INFERIOR - DIVIDIDO) ---
                            col_controles, col_3d = st.columns([1.2, 1])

                            with col_controles:
                                st.subheader("⚙️ Configuração da Modelagem")
                                mutacao_selecionada = st.selectbox("Selecione uma variante para análise 3D:", mutacoes_filtradas['Variante'])

                                linha_selecionada = mutacoes_filtradas[mutacoes_filtradas['Variante'] == mutacao_selecionada]
                                posicao_mutacao = int(linha_selecionada['Posicao'].values[0])

                                if len(mutacao_selecionada) >= 8:
                                    aa_original_3l = mutacao_selecionada[2:5]
                                    aa_mutado_3l = mutacao_selecionada[-3:]
                                    aa_orig_1l = MAPA_AMINOACIDOS.get(aa_original_3l, '?')
                                    aa_mut_1l = MAPA_AMINOACIDOS.get(aa_mutado_3l, '?')

                                    if posicao_mutacao <= len(sequencia_selvagem):
                                        aa_uniprot = sequencia_selvagem[posicao_mutacao - 1]
                                        
                                        if aa_uniprot != aa_orig_1l:
                                            st.warning(f"⚠️ Alerta de Isoforma: A variante espera '{aa_original_3l}' na posição {posicao_mutacao}, mas o UniProt possui '{aa_uniprot}'.")
                                        
                                        sequencia_mutada = sequencia_selvagem[:posicao_mutacao - 1] + aa_mut_1l + sequencia_selvagem[posicao_mutacao:]

                                        if st.button("🚀 Iniciar Modelagem por Homologia (SWISS-MODEL)", use_container_width=True):
                                            if not token_swiss:
                                                st.error("Credencial SWISS-MODEL ausente. Insira o Token no menu lateral.")
                                            else:
                                                with st.spinner("Enviando requisição (POST) para o servidor..."):
                                                    token_limpo = "".join(c for c in token_swiss if c.isalnum() or c == '-')
                                                    headers = {
                                                        "Authorization": f"Token {token_limpo}",
                                                        "Content-Type": "application/json"
                                                    }
                                                    payload = {
                                                        "target_sequences": [sequencia_mutada],
                                                        "project_title": f"GenoStruct_{gene_buscado}_{mutacao_selecionada}"
                                                    }

                                                    try:
                                                        url_swiss = "https://swissmodel.expasy.org/automodel"
                                                        resposta_swiss = requests.post(url_swiss, headers=headers, json=payload, timeout=15)

                                                        if resposta_swiss.status_code in [200, 201, 202]:
                                                            dados_projeto = resposta_swiss.json()
                                                            id_projeto = dados_projeto.get('project_id', 'Desconhecido')
                                                            st.success(f"Submissão autorizada! Job ID: {id_projeto}")
                                                            st.markdown(f"🔗 [Acompanhar renderização no dashboard oficial](https://swissmodel.expasy.org/interactive/)")
                                                        else:
                                                            st.error(f"Erro do Servidor (Código {resposta_swiss.status_code}): Requisição negada.")
                                                    except requests.exceptions.RequestException:
                                                        st.error("Falha de rede ao tentar conectar com a Suíça.")
                                    else:
                                        st.error("Erro: A posição da mutação excede o tamanho da proteína devolvida pelo UniProt.")
                                else:
                                    st.error("Erro de formatação na variante (esperado padrão p.XxxYYYzzz).")
                            
                            with col_3d:
                                st.subheader("🔬 Renderização Estrutural")
                                if posicao_mutacao:
                                    with st.spinner("Baixando coordenadas atômicas PDB..."):
                                        pdb_texto = buscar_pdb_alphafold(uniprot_id)

                                    if pdb_texto:
                                        st.caption(f"Visualizando modelo selvagem. Resíduo Mutado ({posicao_mutacao}) em vermelho.")
                                        view = py3Dmol.view(width=450, height=450)
                                        view.addModel(pdb_texto, 'pdb')
                                        view.setStyle({'cartoon': {'color': 'lightgray'}})
                                        view.setStyle({'resi': str(posicao_mutacao)}, {'stick': {'colorscheme': 'redCarbon', 'radius': 0.3}})
                                        view.zoomTo()
                                        showmol(view, height=450, width=450)
                                    else:
                                        st.warning("Modelo estrutural primário indisponível no repositório AlphaFold.")
                        else:
                            st.info("Nenhuma variante registrada para este alvo no banco de dados local.")
                    else:
                        st.error("⚠️ Alvo genômico não reconhecido nos repositórios oficiais.")
                except requests.exceptions.RequestException:
                    st.error("⚠️ Falha de comunicação com os bancos de dados. Verifique sua conexão de rede.")
else:
    st.info("Arquivo de banco de dados não encontrado na pasta 'dados'.")

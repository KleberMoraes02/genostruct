import streamlit as st
import requests
import pandas as pd
import py3Dmol
from stmol import showmol
import os
import plotly.express as px

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
        caminho = os.path.join("dados", "banco_teste.csv.gz")
        df = pd.read_csv(caminho)
        
        df = df.replace('Sem registro clÃ­nico', 'Sem registro clínico')

        if all(col in df.columns for col in ['Chromosome', 'Position', 'Ref', 'Alt']):
            df['Link_gnomAD'] = "https://gnomad.broadinstitute.org/variant/" + \
                                df['Chromosome'].astype(str) + "-" + \
                                df['Position'].astype(str) + "-" + \
                                df['Ref'].astype(str) + "-" + \
                                df['Alt'].astype(str) + "?dataset=gnomad_r4"
        return df
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
    st.caption("v1.2.1 (Dashboard Integrado)")
    st.markdown("---")
    st.markdown("**⚙️ Credenciais de Modelagem**")
    token_swiss = st.text_input("SWISS-MODEL API Token:", type="password")

# --- CORPO PRINCIPAL DO SITE ---
st.title("GenoStruct: Integração Genômica e Estrutural")

df_mutacoes = carregar_banco_mutacoes()

if not df_mutacoes.empty:
    
    # ==========================================
    # SEÇÃO 1: PESQUISA DO GENE E 3D (TOPO)
    # ==========================================
    st.header("🧬 Análise por Gene e Estrutura 3D")
    st.markdown("Busque um gene para correlacionar variantes clínicas com predições estruturais.")
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

                        st.subheader("📊 Perfil Mutacional Clínico")
                        mutacoes_filtradas = df_mutacoes[df_mutacoes['Gene'].str.upper() == gene_buscado]

                        with st.expander("📖 Dicionário de Colunas (Clique para expandir)"):
                            st.markdown("""
                            Este guia explica os dados genômicos, clínicos e preditivos exibidos na tabela abaixo.
                            
                            **Identificação e Genética:**
                            - **Gene:** Símbolo oficial aprovado pelo comitê HGNC.
                            - **HGNC_ID:** Identificador numérico único do gene.
                            - **UniProt_Accession:** Código do gene no banco mundial de proteínas (UniProt).
                            - **rsID:** Identificador universal da variante no banco dbSNP.
                            
                            **Nomenclatura (Padrão HGVS):**
                            - **HGVS_c:** Alteração no nível do DNA (sequência codificante).
                            - **HGVS_p / Variante:** Alteração no nível da Proteína.
                            - **Posicao:** A posição numérica exata do aminoácido afetado na cadeia proteica.
                            
                            **Predição e Relevância Clínica:**
                            - **AlphaMissense_Class:** Inteligência artificial do Google DeepMind que prevê se a mutação é Patogênica, Benigna ou Ambígua.
                            - **ClinVar:** Registro de relevância clínica real observada em pacientes médicos.
                            - **Link_gnomAD:** Link direto para consultar a frequência da mutação na população mundial (banco gnomAD).
                            """)

                        if not mutacoes_filtradas.empty:
                            todas_as_colunas = mutacoes_filtradas.columns.tolist()
                            colunas_padrao = ['Gene', 'Variante', 'Posicao', 'Troca', 'AlphaMissense_Class', 'ClinVar', 'Link_gnomAD']
                            colunas_padrao = [c for c in colunas_padrao if c in todas_as_colunas]

                            colunas_selecionadas = st.multiselect(
                                "Selecione as colunas que deseja visualizar na tabela:",
                                options=todas_as_colunas,
                                default=colunas_padrao
                            )

                            df_exibicao = mutacoes_filtradas[colunas_selecionadas]
                            
                            configuracao_colunas = {}
                            if 'Link_gnomAD' in df_exibicao.columns:
                                configuracao_colunas['Link_gnomAD'] = st.column_config.LinkColumn("🔗 Frequência (gnomAD)", display_text="Ver Variante")

                            st.dataframe(df_exibicao, use_container_width=True, hide_index=True, column_config=configuracao_colunas)
                            st.markdown("---")
                            
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
                                                    headers = {"Authorization": f"Token {token_limpo}", "Content-Type": "application/json"}
                                                    payload = {"target_sequences": [sequencia_mutada], "project_title": f"GenoStruct_{gene_buscado}_{mutacao_selecionada}"}

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

    # ==========================================
    # SEÇÃO 2: GRÁFICOS DO BANCO (BAIXO)
    # ==========================================
    st.markdown("<br><br>", unsafe_allow_html=True) # Dá um espaço em branco extra
    st.markdown("---")
    st.header("📈 Visão Geral do Banco de Dados")
    st.markdown(f"**Total de variantes carregadas no banco:** {len(df_mutacoes)}")
    
    col_graf_1, col_graf_2 = st.columns(2)
    
    with col_graf_1:
        if 'Gene' in df_mutacoes.columns:
            st.subheader("Top 10 Genes com mais Mutações")
            top_genes = df_mutacoes['Gene'].value_counts().head(10).reset_index()
            top_genes.columns = ['Gene', 'Quantidade']
            
            fig_bar = px.bar(top_genes, x='Gene', y='Quantidade', color='Quantidade', color_continuous_scale='Blues')
            st.plotly_chart(fig_bar, use_container_width=True)
            
    with col_graf_2:
        if 'AlphaMissense_Class' in df_mutacoes.columns:
            st.subheader("Predição AlphaMissense Global")
            am_counts = df_mutacoes['AlphaMissense_Class'].value_counts().reset_index()
            am_counts.columns = ['Classificação', 'Total']
            
            cores_am = {'Pathogenic':'#ff4b4b', 'Benign':'#1f77b4', 'Ambiguous':'#ffc107', 'Não avaliado':'#d3d3d3'}
            fig_pie = px.pie(am_counts, names='Classificação', values='Total', hole=0.4, color='Classificação', color_discrete_map=cores_am)
            st.plotly_chart(fig_pie, use_container_width=True)

else:
    st.info("Arquivo de banco de dados não encontrado na pasta 'dados'.")

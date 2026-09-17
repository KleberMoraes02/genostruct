import streamlit as st
import requests
import pandas as pd
import py3Dmol
from stmol import showmol
import os
import plotly.express as px
import time

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
    st.caption("v1.5.0 (IA Biológica e Automação)")
    st.markdown("---")
    st.markdown("**⚙️ Credenciais de Modelagem**")
    token_swiss = st.text_input("SWISS-MODEL API Token:", type="password")

# --- CORPO PRINCIPAL DO SITE ---
st.title("GenoStruct: Integração Genômica e Estrutural")

df_mutacoes = carregar_banco_mutacoes()

if not df_mutacoes.empty:
    
    gene_buscado_raw = st.text_input("Nome do Gene alvo (ex: ABCD2):")
    gene_buscado = gene_buscado_raw.strip().upper()

    if gene_buscado:
        if len(gene_buscado) < 2:
            st.warning("⚠️ Digite um nome de gene válido.")
        else:
            with st.spinner("Sincronizando com UniProtKB..."):
                try: 
                    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_buscado} AND organism_id:9606 AND reviewed:true&format=json&size=1"
                    resposta = requests.get(url, timeout=10) 

                    if resposta.status_code == 200 and len(resposta.json()['results']) > 0:
                        proteina = resposta.json()['results'][0]
                        uniprot_id = proteina['primaryAccession']
                        sequencia_selvagem = proteina['sequence']['value']
                        
                        # NOVO: EXTRAINDO MAPA DE CARACTERÍSTICAS DA PROTEÍNA
                        mapa_funcional = proteina.get('features', [])

                        st.success(f"Conexão estabelecida! Alvo: **{gene_buscado}** (Accession: {uniprot_id})")

                        mutacoes_filtradas = df_mutacoes[df_mutacoes['Gene'].str.upper() == gene_buscado]
                        st.subheader("📊 Perfil Mutacional Clínico")

                        with st.expander("📖 Dicionário de Colunas"):
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

                            colunas_selecionadas = st.multiselect("Selecione as colunas:", options=todas_as_colunas, default=colunas_padrao)
                            df_exibicao = mutacoes_filtradas[colunas_selecionadas]
                            
                            configuracao_colunas = {}
                            if 'Link_gnomAD' in df_exibicao.columns:
                                configuracao_colunas['Link_gnomAD'] = st.column_config.LinkColumn("🔗 Frequência", display_text="Ver Variante")
                            st.dataframe(df_exibicao, use_container_width=True, hide_index=True, column_config=configuracao_colunas)
                            
                            csv_data = df_exibicao.to_csv(sep='\t', index=False).encode('utf-8')
                            st.download_button("📥 Exportar Tabela Curada (.tsv)", data=csv_data, file_name=f"GenoStruct_{gene_buscado}_variantes.tsv", mime="text/tsv")

                            st.markdown("---")
                            
                            # GRÁFICO DE DISPERSÃO
                            st.subheader("📍 Paisagem Mutacional (Mapeamento de Patogenicidade)")
                            if 'AlphaMissense_Class' in mutacoes_filtradas.columns and 'Posicao' in mutacoes_filtradas.columns:
                                cores_am = {'likely_pathogenic': '#ff4b4b', 'likely_benign': '#1f77b4', 'ambiguous': '#ffc107', 'Não disponível': '#d3d3d3'}
                                fig_scatter = px.scatter(mutacoes_filtradas, x='Posicao', y='AlphaMissense_Class', color='AlphaMissense_Class', color_discrete_map=cores_am, hover_data=['Variante', 'Troca', 'ClinVar'])
                                fig_scatter.update_traces(marker=dict(size=10, line=dict(width=1, color='DarkSlateGrey')))
                                st.plotly_chart(fig_scatter, use_container_width=True)
                            
                            st.markdown("---")
                            
                            # CONTROLES E ESTRUTURA 3D
                            col_controles, col_3d = st.columns([1.2, 1])

                            with col_controles:
                                st.subheader("⚙️ Configuração da Modelagem")
                                mutacao_selecionada = st.selectbox("Selecione uma variante para análise 3D:", mutacoes_filtradas['Variante'])
                                
                                # Limpa o modelo 3D mutado se o usuário trocar de mutação
                                seletor_mudou = 'ultima_mutacao' not in st.session_state or st.session_state['ultima_mutacao'] != mutacao_selecionada
                                if seletor_mudou:
                                    st.session_state['pdb_mutante'] = None
                                    st.session_state['ultima_mutacao'] = mutacao_selecionada

                                linha_selecionada = mutacoes_filtradas[mutacoes_filtradas['Variante'] == mutacao_selecionada]
                                posicao_mutacao = int(linha_selecionada['Posicao'].values[0])

                                if len(mutacao_selecionada) >= 8:
                                    aa_original_3l = mutacao_selecionada[2:5]
                                    aa_mutado_3l = mutacao_selecionada[-3:]
                                    aa_orig_1l = MAPA_AMINOACIDOS.get(aa_original_3l, '?')
                                    aa_mut_1l = MAPA_AMINOACIDOS.get(aa_mutado_3l, '?')

                                    # INTELIGÊNCIA BIOLÓGICA (UniProt Features)
                                    is_peptideo_sinal = False
                                    for feature in mapa_funcional:
                                        tipo_alvo = feature.get('type', '')
                                        try:
                                            inicio = int(feature['location']['start']['value'])
                                            fim = int(feature['location']['end']['value'])
                                            
                                            # Verifica se a mutação caiu no meio dessa característica
                                            if inicio <= posicao_mutacao <= fim:
                                                if tipo_alvo == 'Signal':
                                                    is_peptideo_sinal = True
                                                    st.info(f"🟡 **Peptídeo Sinal:** A mutação cai na região de clivagem (pos {inicio}-{fim}). É normal haver divergência de isoforma aqui.")
                                                elif tipo_alvo == 'Active site':
                                                    st.error(f"🔴 **SÍTIO ATIVO DESTRUÍDO:** A mutação afeta o sítio catalítico crítico da proteína na posição {posicao_mutacao}!")
                                                elif tipo_alvo == 'Disulfide bond':
                                                    st.warning(f"🟠 **Ponte Dissulfeto Ameaçada:** Esta mutação pode quebrar uma ligação estrutural importante.")
                                        except:
                                            pass

                                    if posicao_mutacao <= len(sequencia_selvagem):
                                        aa_uniprot = sequencia_selvagem[posicao_mutacao - 1]
                                        
                                        if aa_uniprot != aa_orig_1l and not is_peptideo_sinal:
                                            st.warning(f"⚠️ Alerta de Isoforma: Variante espera '{aa_original_3l}', mas o UniProt possui '{aa_uniprot}'.")
                                        
                                        sequencia_mutada = sequencia_selvagem[:posicao_mutacao - 1] + aa_mut_1l + sequencia_selvagem[posicao_mutacao:]

                                        # AUTOMAÇÃO DO SWISS-MODEL
                                        if st.button("🚀 Iniciar Automação SWISS-MODEL", use_container_width=True):
                                            if not token_swiss:
                                                st.error("Credencial SWISS-MODEL ausente no menu lateral.")
                                            else:
                                                token_limpo = "".join(c for c in token_swiss if c.isalnum() or c == '-')
                                                headers = {"Authorization": f"Token {token_limpo}", "Content-Type": "application/json"}
                                                payload = {"target_sequences": [sequencia_mutada], "project_title": f"GenoStruct_{mutacao_selecionada}"}

                                                try:
                                                    url_post = "https://swissmodel.expasy.org/automodel"
                                                    resposta_swiss = requests.post(url_post, headers=headers, json=payload, timeout=15)

                                                    if resposta_swiss.status_code in [200, 201, 202]:
                                                        id_projeto = resposta_swiss.json().get('project_id')
                                                        
                                                        # Interface de carregamento dinâmico
                                                        painel_status = st.empty()
                                                        barra_progresso = st.progress(0)
                                                        
                                                        url_status = f"https://swissmodel.expasy.org/project/{id_projeto}/models/summary"
                                                        status = "PENDING"
                                                        tentativas = 0
                                                        
                                                        # LOOP DE ESPERA (Consulta a cada 10s. Máx 3 min)
                                                        while status in ["PENDING", "RUNNING", "QUEUED"] and tentativas < 18:
                                                            time.sleep(10)
                                                            tentativas += 1
                                                            try:
                                                                res_status = requests.get(url_status, headers={"Authorization": f"Token {token_limpo}"}, timeout=10)
                                                                if res_status.status_code == 200:
                                                                    dados_status = res_status.json()
                                                                    status = dados_status.get("status", "UNKNOWN")
                                                                    
                                                                    if status == "QUEUED":
                                                                        painel_status.info(f"⏳ Na fila do servidor suíço... (Tentativa {tentativas}/18)")
                                                                        barra_progresso.progress(20)
                                                                    elif status == "RUNNING":
                                                                        painel_status.warning(f"⚙️ Construindo conformação 3D... (Tentativa {tentativas}/18)")
                                                                        barra_progresso.progress(60)
                                                                    elif status == "COMPLETED":
                                                                        painel_status.success("✅ Modelagem Concluída!")
                                                                        barra_progresso.progress(100)
                                                                        
                                                                        # BAIXANDO O ARQUIVO PDB MUTADO AUTOMATICAMENTE
                                                                        url_pdb = f"https://swissmodel.expasy.org/project/{id_projeto}/models/01.pdb"
                                                                        res_pdb = requests.get(url_pdb, headers={"Authorization": f"Token {token_limpo}"}, timeout=10)
                                                                        if res_pdb.status_code == 200:
                                                                            st.session_state['pdb_mutante'] = res_pdb.text
                                                                    elif status in ["FAILED", "CANCELLED"]:
                                                                        painel_status.error("❌ A modelagem falhou no servidor do SWISS-MODEL.")
                                                                        break
                                                            except:
                                                                painel_status.error("Falha ao checar status. Tentando novamente...")
                                                                
                                                        if tentativas >= 18:
                                                            painel_status.error("Tempo limite de 3 minutos excedido. O servidor está lotado.")
                                                    else:
                                                        st.error(f"Erro do Servidor (Código {resposta_swiss.status_code})")
                                                except requests.exceptions.RequestException:
                                                    st.error("Falha de rede ao conectar com a Suíça.")
                                    else:
                                        st.error("Posição excede o tamanho da proteína.")
                                
                                # SE O DOWNLOAD DEU CERTO, MOSTRA O MUTANTE AQUI NA ESQUERDA
                                if st.session_state.get('pdb_mutante'):
                                    st.markdown("---")
                                    st.subheader("🔴 Estrutura Mutada (SWISS-MODEL)")
                                    view_mut = py3Dmol.view(width=450, height=450)
                                    view_mut.addModel(st.session_state['pdb_mutante'], 'pdb')
                                    view_mut.setStyle({'cartoon': {'color': 'lightblue'}})
                                    view_mut.setStyle({'resi': str(posicao_mutacao)}, {'stick': {'colorscheme': 'redCarbon', 'radius': 0.3}})
                                    view_mut.zoomTo()
                                    showmol(view_mut, height=450, width=450)
                                    
                                    st.download_button(
                                        label="📥 Exportar Mutante Modelado (.pdb)",
                                        data=st.session_state['pdb_mutante'],
                                        file_name=f"GenoStruct_{gene_buscado}_{mutacao_selecionada}.pdb",
                                        mime="chemical/x-pdb",
                                    )

                            with col_3d:
                                st.subheader("🔵 Estrutura Selvagem (AlphaFold)")
                                if posicao_mutacao:
                                    with st.spinner("Baixando coordenadas PDB..."):
                                        pdb_texto = buscar_pdb_alphafold(uniprot_id)
                                    if pdb_texto:
                                        view = py3Dmol.view(width=450, height=450)
                                        view.addModel(pdb_texto, 'pdb')
                                        view.setStyle({'cartoon': {'color': 'lightgray'}})
                                        view.setStyle({'resi': str(posicao_mutacao)}, {'stick': {'colorscheme': 'redCarbon', 'radius': 0.3}})
                                        view.zoomTo()
                                        showmol(view, height=450, width=450)
                                        
                                        st.download_button(
                                            label="📥 Exportar Selvagem (.pdb)",
                                            data=pdb_texto,
                                            file_name=f"AlphaFold_{uniprot_id}.pdb",
                                            mime="chemical/x-pdb",
                                        )

                        # AQUI ESTÃO AS LINHAS QUE TINHAM SIDO CORTADAS:
                        else:
                            st.info("Nenhuma variante registrada no banco local para este gene.")
                    else:
                        st.error("⚠️ Alvo genômico não reconhecido nos repositórios oficiais.")
                except requests.exceptions.RequestException:
                    st.error("⚠️ Falha de comunicação de rede ao tentar contatar o UniProt.")
else:
    st.info("Arquivo de banco de dados não encontrado na pasta 'dados'.")

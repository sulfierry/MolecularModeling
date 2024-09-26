# Análise de Interações Proteicas

Este conjunto de scripts foi desenvolvido para analisar interações entre cadeias de proteínas a partir de simulações de dinâmica molecular. O processo envolve a extração de dados de trajetória, processamento de interações e visualização dos resultados.

## Conteúdo

1. `near_atom_from_chain.sh`: Script principal em Bash que coordena todo o processo de análise.
2. `process_interaction.py`: Script Python gerado dinamicamente para processar e visualizar as interações.

## Requisitos

- Bash
- Python 3.x
- VMD (Visual Molecular Dynamics)
- Bibliotecas Python: pandas, matplotlib, seaborn

## Instalação

1. Clone este repositório ou baixe os scripts para seu computador.
2. Certifique-se de que o VMD está instalado e acessível pelo terminal.
3. Instale as bibliotecas Python necessárias:
   ```
   pip install pandas matplotlib seaborn
   ```

## Uso

1. Coloque seus arquivos de trajetória (`.dcd`) e topologia (`.prmtop`) na mesma pasta que os scripts.

2. Abra o script `near_atom_from_chain.sh` e ajuste as seguintes variáveis conforme necessário:
   - `BASE_NAME`: Nome base dos seus arquivos de simulação
   - `TOPOLOGY`: Nome do arquivo de topologia
   - `COORDINATES`: Nome do arquivo de trajetória
   - `TOTAL_FRAMES`: Número total de frames na sua trajetória
   - `CHAIN_A_RESIDUES` e `CHAIN_B_RESIDUES`: Definição dos resíduos para cada cadeia
   - `CUTOFF_DISTANCE`: Distância de corte para considerar uma interação
   - `TRESHOLD_PREVALENCE_INTERACTION`: Limiar de prevalência para filtrar interações

3. Execute o script Bash:
   ```
   bash near_atom_from_chain.sh
   ```

4. O script irá processar os dados e gerar vários arquivos de saída, incluindo:
   - `detailed_interactions.tsv`: Dados brutos de todas as interações
   - `aminoacid_interactions_prevalence.csv`: Dados de prevalência das interações
   - `prev_interaction_top_X.png`: Gráfico das X interações mais prevalentes
   - `prevalence.png`: Gráfico de barras mostrando a prevalência de todas as interações

## Notas Importantes

- Certifique-se de que a pasta `vmd_frames` está no diretório atual antes de executar o script. Esta pasta é criada após a execução do script `near_atom_from_chain.sh`.
- O script processa os frames em blocos para otimizar o uso de memória. Você pode ajustar o tamanho do bloco modificando a variável `SLICE` no script Bash.
- A análise pode levar um tempo considerável dependendo do tamanho da sua trajetória e do número de interações.
- Os resultados são salvos em arquivos CSV e imagens PNG para fácil visualização e análise posterior.

## Solução de Problemas

- Se encontrar erros relacionados ao VMD, verifique se ele está corretamente instalado e acessível pelo terminal.
- Para problemas com bibliotecas Python, certifique-se de que todas as dependências estão instaladas e atualizadas.
- Se o script parar inesperadamente, verifique se há espaço suficiente em disco e se você tem permissões de escrita no diretório atual.

## Contribuições

Contribuições para melhorar estes scripts são bem-vindas. Por favor, abra uma issue ou envie um pull request com suas sugestões ou melhorias.

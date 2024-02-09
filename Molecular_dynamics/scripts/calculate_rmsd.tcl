# Nome do script: calcular_rmsd_alinhado_corrigido.tcl

# Definindo os arquivos de entrada
set topologia "5cc8_wr_1.prmtop"
set trajetoria "5cc8_wr_1.dcd"

# Carregando a topologia e a trajetória
mol new $topologia type parm7
mol addfile $trajetoria type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

# Selecionando os átomos do backbone e seus átomos pesados (C, CA, O, N)
set selecao [atomselect top "backbone and (name C or name CA or name O or name N)"]

# Inicializa o arquivo de saída para os resultados do RMSD
set arquivo_saida "rmsd_alinhado_resultados_corrigido.txt"
set fileId [open $arquivo_saida "w"]

# Alinhando a trajetória e calculando o RMSD
set totalFrames [molinfo top get numframes]
set refCoords ""
for {set frame 0} {$frame < $totalFrames} {incr frame} {
    $selecao frame $frame
    $selecao update
    if {$frame == 0} {
        # Guarda as coordenadas do primeiro frame para usar como referência
        set refCoords [$selecao get {x y z}]
    } else {
        # Alinha cada frame com o primeiro frame
        set moveMatrix [measure fit $selecao $refCoords]
        $selecao move $moveMatrix
    }
    # Calcula o RMSD do frame alinhado em relação ao primeiro frame
    set rmsd [measure rmsd $selecao $refCoords]
    puts $fileId "$frame $rmsd"
}

# Fechando o arquivo de saída
close $fileId

# Finalizando o script
exit

# Nome do script: calcular_rmsd_alinhado_corrigido.tcl

# Definindo os arquivos de entrada
set topologia "5cc8_wr_1.prmtop"
set trajetoria "5cc8_wr_1.dcd"

# Carregando a topologia e a trajetória
mol new $topologia type parm7
mol addfile $trajetoria type dcd first 0 last -1 step 1 waitfor all

# Selecionando os átomos do backbone e seus respectivos átomos pesados (C, CA, O, N)
set selecao "backbone and (name C or name CA or name O or name N)"

# Inicializa o arquivo de saída para os resultados do RMSD
set arquivo_saida "rmsd_alinhado_resultados_corrigido.txt"
set fileId [open $arquivo_saida "w"]

# Prepara a estrutura de referência (primeiro frame)
set ref [atomselect top "$selecao" frame 0]
set refCoords [$ref get {x y z}]

# Alinhando a trajetória e calculando o RMSD
set totalFrames [molinfo top get numframes]
for {set frame 0} {$frame < $totalFrames} {incr frame} {
    set sel [atomselect top "$selecao" frame $frame]
    if {$frame > 0} {
        # Alinha cada frame com o primeiro frame
        set moveMatrix [measure fit $sel $refCoords]
        $sel move $moveMatrix
    }
    # Calcula o RMSD do frame alinhado em relação ao primeiro frame
    set rmsd [measure rmsd $sel $ref]
    puts $fileId "$frame $rmsd"
    $sel delete
}

# Fechamento da seleção de referência e do arquivo de saída
$ref delete
close $fileId

# Finalizando o script
exit

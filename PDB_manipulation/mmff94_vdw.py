#                Copyright (c) Merck and Co., Inc., 1994
#                         All Rights Reserved

# Transcribed by: Leon Sulfierry GitHub:https://github.com/sulfierry

# 'Valor_escalar' : Esta é a chave do dicionário, que provavelmente indica um tipo de átomo específico ou uma classificação MMFF94 para um átomo.
# 'alpha-i' :Representa o parâmetro ϵ no potencial de Lennard-Jones. É a profundidade do poço de potencial e descreve a força de atração entre os átomos.
# 'N-i' :  Equivale ao parâmetro σ no potencial de Lennard-Jones. É a distância na qual o potencial entre dois átomos é zero, basicamente referindo-se ao tamanho efetivo do átomo.
# 'A-i' e 'G-i': Estes estão associados ao termo de energia de ângulo de ligação no MMFF94.
# 'A-i': Refere-se ao ângulo de ligação de equilíbrio em graus.
# 'G-i': É a constante de força associada ao ângulo de ligação.
# 'DA': Em muitos campos de força, "DA" refere-se a um átomo doador/aceitador. No contexto do MMFF94, pode indicar se o átomo em questão é um doador ou aceitador de ligações de hidrogênio ou pode ter algum outro significado específico.
# 'Symb': Representa um símbolo ou abreviação para o tipo de átomo ou grupo funcional.
# Origin': Indica a origem ou fonte dos parâmetros. "E94" provavelmente se refere a uma especificação ou literatura específica de onde os parâmetros foram derivados ou adaptados para o MMFF94.

mmff94_vdw = {
    '1': {'alpha-i': '1.050', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CR', 'Origin': 'E94'},
    '2': {'alpha-i': '1.350', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'C=C', 'Origin': 'E94'},
    '3': {'alpha-i': '1.100', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'C=O', 'Origin': 'E94'},
    '4': {'alpha-i': '1.300', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CSP', 'Origin': 'E94'},
    '5': {'alpha-i': '0.250', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': '-', 'Symb': 'HC', 'Origin': 'C94'},
    '6': {'alpha-i': '0.70', 'N-i': '3.150', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'OR', 'Origin': 'C94'},
    '7': {'alpha-i': '0.65', 'N-i': '3.150', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'O=C', 'Origin': 'C94'},
    '8': {'alpha-i': '1.15', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NR', 'Origin': 'C94'},
    '9': {'alpha-i': '0.90', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'N=C', 'Origin': 'C94'},
    '10': {'alpha-i': '1.000', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NC=O', 'Origin': 'E94'},
    '11': {'alpha-i': '0.35', 'N-i': '3.480', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'F', 'Origin': 'C94'},
    '12': {'alpha-i': '2.300', 'N-i': '5.100', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': 'CL', 'Origin': 'E94'},
    '13': {'alpha-i': '3.400', 'N-i': '6.000', 'A-i': '3.190', 'G-i': '1.359', 'DA': 'A', 'Symb': 'BR', 'Origin': 'E94'},
    '14': {'alpha-i': '5.500', 'N-i': '6.950', 'A-i': '3.080', 'G-i': '1.404', 'DA': 'A', 'Symb': 'I', 'Origin': 'E94'},
    '15': {'alpha-i': '3.00', 'N-i': '4.800', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': 'S', 'Origin': 'C94'},
    '16': {'alpha-i': '3.900', 'N-i': '4.800', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': 'S=C', 'Origin': 'E94'},
    '17': {'alpha-i': '2.700', 'N-i': '4.800', 'A-i': '3.320', 'G-i': '1.345', 'DA': '-', 'Symb': 'SO', 'Origin': 'E94'},
    '18': {'alpha-i': '2.100', 'N-i': '4.800', 'A-i': '3.320', 'G-i': '1.345', 'DA': '-', 'Symb': 'SO2', 'Origin': 'E94'},
    '19': {'alpha-i': '4.500', 'N-i': '4.200', 'A-i': '3.320', 'G-i': '1.345', 'DA': '-', 'Symb': 'SI', 'Origin': 'E94'},
    '20': {'alpha-i': '1.050', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CR3R', 'Origin': 'E94'},
    '21': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HOR', 'Origin': 'C94'},
    '22': {'alpha-i': '1.100', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CR3R', 'Origin': 'E94'},
    '23': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HNR', 'Origin': 'C94'},
    '24': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HOCO', 'Origin': 'C94'},
    '25': {'alpha-i': '1.600', 'N-i': '4.500', 'A-i': '3.320', 'G-i': '1.345', 'DA': '-', 'Symb': 'PO4', 'Origin': 'E94'},
    '26': {'alpha-i': '3.600', 'N-i': '4.500', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': 'P', 'Origin': 'E94'},
    '27': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HN=C', 'Origin': 'C94'},
    '28': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HNCO', 'Origin': 'C94'},
    '29': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HOCC', 'Origin': 'C94'},
    '30': {'alpha-i': '1.350', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CE4R', 'Origin': 'E94'},
    '31': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HOH', 'Origin': 'C94'},
    '32': {'alpha-i': '0.75', 'N-i': '3.150', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'O2CM', 'Origin': 'C94'},
    '33': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HOS', 'Origin': 'C94'},
    '34': {'alpha-i': '1.00', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'NR+', 'Origin': 'C94'},
    '35': {'alpha-i': '1.50', 'N-i': '3.150', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'OM', 'Origin': 'X94'},
    '36': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HNR+', 'Origin': 'C94'},
    '37': {'alpha-i': '1.350', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CB', 'Origin': 'E94'},
    '38': {'alpha-i': '0.85', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NPYD', 'Origin': 'C94'},
    '39': {'alpha-i': '1.10', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'NPYL', 'Origin': 'C94'},
    '40': {'alpha-i': '1.00', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NC=C', 'Origin': 'E94'},
    '41': {'alpha-i': '1.100', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CO2M', 'Origin': 'C94'},
    '42': {'alpha-i': '1.000', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NSP', 'Origin': 'E94'},
    '43': {'alpha-i': '1.000', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NSO2', 'Origin': 'E94'},
    '44': {'alpha-i': '3.00', 'N-i': '4.800', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': 'STHI', 'Origin': 'C94'},
    '45': {'alpha-i': '1.150', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'NO2', 'Origin': 'E94'},
    '46': {'alpha-i': '1.300', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'N=O', 'Origin': 'E94'},
    '47': {'alpha-i': '1.000', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NAZT', 'Origin': 'X94'},
    '48': {'alpha-i': '1.200', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NSO', 'Origin': 'X94'},
    '49': {'alpha-i': '1.00', 'N-i': '3.150', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'O+', 'Origin': 'X94'},
    '50': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HO+', 'Origin': 'C94'},
    '51': {'alpha-i': '0.400', 'N-i': '3.150', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'O=+', 'Origin': 'E94'},
    '52': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HO=+', 'Origin': 'C94'},
    '53': {'alpha-i': '1.000', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': '=N=', 'Origin': 'X94'},
    '54': {'alpha-i': '1.30', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'N+=C', 'Origin': 'C94'},
    '55': {'alpha-i': '0.80', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'NCN+', 'Origin': 'E94'},
    '56': {'alpha-i': '0.80', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'NGD+', 'Origin': 'E94'},
    '57': {'alpha-i': '1.000', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CNN+', 'Origin': 'E94'},
    '58': {'alpha-i': '0.80', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'NPD+', 'Origin': 'E94'},
    '59': {'alpha-i': '0.65', 'N-i': '3.150', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'OFUR', 'Origin': 'C94'},
    '60': {'alpha-i': '1.800', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'C%-', 'Origin': 'E94'},
    '61': {'alpha-i': '0.800', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NR%', 'Origin': 'E94'},
    '62': {'alpha-i': '1.300', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NM', 'Origin': 'X94'},
    '63': {'alpha-i': '1.350', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'C5A', 'Origin': 'E94'},
    '64': {'alpha-i': '1.350', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'C5B', 'Origin': 'E94'},
    '65': {'alpha-i': '1.000', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'N5A', 'Origin': 'E94'},
    '66': {'alpha-i': '0.75', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'N5B', 'Origin': 'C94'},
    '67': {'alpha-i': '0.950', 'N-i': '2.82', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'N2OX', 'Origin': 'X94'},
    '68': {'alpha-i': '0.90', 'N-i': '2.82', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'N3OX', 'Origin': 'C94'},
    '69': {'alpha-i': '0.950', 'N-i': '2.82', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NPOX', 'Origin': 'C94'},
    '70': {'alpha-i': '0.87', 'N-i': '3.150', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'OH2', 'Origin': 'C94'},
    '71': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HS', 'Origin': 'C94'},
    '72': {'alpha-i': '4.000', 'N-i': '4.800', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': 'SM', 'Origin': 'X94'},
    '73': {'alpha-i': '3.000', 'N-i': '4.800', 'A-i': '3.320', 'G-i': '1.345', 'DA': '-', 'Symb': 'SMO2', 'Origin': 'X94'},
    '74': {'alpha-i': '3.000', 'N-i': '4.800', 'A-i': '3.320', 'G-i': '1.345', 'DA': '-', 'Symb': '=S=O', 'Origin': 'X94'},
    '75': {'alpha-i': '4.000', 'N-i': '4.500', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': '-P=C', 'Origin': 'X94'},
    '76': {'alpha-i': '1.200', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'N5M', 'Origin': 'X94'},
    '77': {'alpha-i': '1.500', 'N-i': '5.100', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': 'CLO4', 'Origin': 'X94'},
    '78': {'alpha-i': '1.350', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'C5', 'Origin': 'X94'},
    '79': {'alpha-i': '1.000', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'N5', 'Origin': 'X94'},
    '80': {'alpha-i': '1.000', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CIM+', 'Origin': 'C94'},
    '81': {'alpha-i': '0.80', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'NIM+', 'Origin': 'C94'},
    '82': {'alpha-i': '0.950', 'N-i': '2.82', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'N5OX', 'Origin': 'X94'},
    '87': {'alpha-i': '0.45', 'N-i': '6.', 'A-i': '4.', 'G-i': '1.4', 'DA': '-', 'Symb': 'FE+2', 'Origin': 'X94'},
    '88': {'alpha-i': '0.55', 'N-i': '6.', 'A-i': '4.', 'G-i': '1.4', 'DA': '-', 'Symb': 'FE+3', 'Origin': 'X94'},
    '89': {'alpha-i': '1.4', 'N-i': '3.48', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'F-', 'Origin': 'X94'},
    '90': {'alpha-i': '4.5', 'N-i': '5.100', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': 'CL-', 'Origin': 'X94'},
    '91': {'alpha-i': '6.0', 'N-i': '6.000', 'A-i': '3.190', 'G-i': '1.359', 'DA': 'A', 'Symb': 'BR-', 'Origin': 'X94'},
    '92': {'alpha-i': '0.15', 'N-i': '2.', 'A-i': '4.', 'G-i': '1.3', 'DA': '-', 'Symb': 'LI+', 'Origin': 'X94'},
    '93': {'alpha-i': '0.4', 'N-i': '3.5', 'A-i': '4.', 'G-i': '1.3', 'DA': '-', 'Symb': 'NA+', 'Origin': 'X94'},
    '94': {'alpha-i': '1.0', 'N-i': '5.', 'A-i': '4.', 'G-i': '1.3', 'DA': '-', 'Symb': 'K+', 'Origin': 'X94'},
    '95': {'alpha-i': '0.43', 'N-i': '6.', 'A-i': '4.', 'G-i': '1.4', 'DA': '-', 'Symb': 'ZN+2', 'Origin': 'X94'},
    '96': {'alpha-i': '0.9', 'N-i': '5.', 'A-i': '4.', 'G-i': '1.4', 'DA': '-', 'Symb': 'CA+2', 'Origin': 'X94'},
    '97': {'alpha-i': '0.35', 'N-i': '6.', 'A-i': '4.', 'G-i': '1.4', 'DA': '-', 'Symb': 'CU+1', 'Origin': 'X94'},
    '98': {'alpha-i': '0.40', 'N-i': '6.', 'A-i': '4.', 'G-i': '1.4', 'DA': '-', 'Symb': 'CU+2', 'Origin': 'X94'},
    '99': {'alpha-i': '0.35', 'N-i': '3.5', 'A-i': '4.', 'G-i': '1.3', 'DA': '-', 'Symb': 'MG+2', 'Origin': 'X94'},
}

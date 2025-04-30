#!MC 1410

# Define the variable names to be used
$!VarSet |var_output| = 'LESIQ_f_spanavg'
$!VarSet |var_input|  = 'LES_IQ_f'

$!AlterData 
  IgnoreDivideByZero = Yes
  Equation = '{|var_output|} = 0'

$! Loop |NumZones|
$! AlterData
   IgnoreDivideByZero = Yes
   Equation = '{|var_output|} = {|var_output|} + {|var_input|}[|Loop|]'

$! EndLoop 

$! AlterData
   IgnoreDivideByZero = Yes
   Equation = '{|var_output|} = {|var_output|} / |NumZones|'

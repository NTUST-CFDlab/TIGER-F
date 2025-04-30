#!MC 1410

# Define the variable names to be used
$!VarSet |var_output| = 'result_zoneavg'
$!VarSet |var_input|  = 'input_var'

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

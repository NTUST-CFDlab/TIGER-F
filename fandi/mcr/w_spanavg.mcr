#!MC 1410

# Define the variable names to be used
$!VarSet |var_output| = 'w_spanavg'
$!VarSet |var_input|  = 'w_timeavg'

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

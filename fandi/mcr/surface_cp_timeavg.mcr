#!MC 1410
$!AlterData 
  IgnoreDivideByZero = Yes
  Equation = '{cp_avg}= 0'
$! Loop |NumZones|
$! AlterData
   IgnoreDivideByZero = Yes
   Equation = '{cp_avg} = {cp_avg} +  {C<sub>p</sub>}[|Loop|]'

$! EndLoop 
$! AlterData
     IgnoreDivideByZero = Yes
	 Equation = '{cp_avg} = {cp_avg} / |NumZones|'
#!MC 1410
$!AlterData 
  IgnoreDivideByZero = Yes
  Equation = '{w_timeavg}= 0'
$! Loop |NumZones|
$! AlterData
   IgnoreDivideByZero = Yes
   Equation = '{w_timeavg} = {w_timeavg} +  {RHO-W}[|Loop|]'

$! EndLoop 
$! AlterData
     IgnoreDivideByZero = Yes
	 Equation = '{w_timeavg} = {w_timeavg} / |NumZones|'
#!MC 1410
$!AlterData 
  IgnoreDivideByZero = Yes
  Equation = '{p_avg}= 0'
$! Loop |NumZones|
$! AlterData
   IgnoreDivideByZero = Yes
   Equation = '{p_avg} = {p_avg} +  {RHO}[|Loop|]'

$! EndLoop 
$! AlterData
     IgnoreDivideByZero = Yes
	 Equation = '{p_avg} = {p_avg} / |NumZones|'

$!AlterData 
  IgnoreDivideByZero = Yes
  Equation = '{w_avg}= 0'
$! Loop |NumZones|
$! AlterData
   IgnoreDivideByZero = Yes
   Equation = '{w_avg} = {w_avg} +  {RHO-W}[|Loop|]'

$! EndLoop 
$! AlterData
     IgnoreDivideByZero = Yes
	 Equation = '{w_avg} = {w_avg} / |NumZones|'
	 
$!AlterData 
  IgnoreDivideByZero = Yes
  Equation = '{v_avg}= 0'
$! Loop |NumZones|
$! AlterData
   IgnoreDivideByZero = Yes
   Equation = '{v_avg} = {v_avg} +  {RHO-V}[|Loop|]'

$! EndLoop 
$! AlterData
     IgnoreDivideByZero = Yes
	 Equation = '{v_avg} = {v_avg} / |NumZones|'
	 
$!AlterData 
  IgnoreDivideByZero = Yes
  Equation = '{u_avg}= 0'
$! Loop |NumZones|
$! AlterData
   IgnoreDivideByZero = Yes
   Equation = '{u_avg} = {u_avg} +  {RHO-U}[|Loop|]'

$! EndLoop 
$! AlterData
     IgnoreDivideByZero = Yes
	 Equation = '{u_avg} = {u_avg} / |NumZones|'
	 
  
$!AlterData 
  IgnoreDivideByZero = Yes
  Equation = '{TKE_avg}= 0'
$! Loop |NumZones|
$! AlterData
   IgnoreDivideByZero = Yes
   Equation = '{TKE_avg} = {TKE_avg} +  {E}[|Loop|]'

$! EndLoop
$! AlterData
     IgnoreDivideByZero = Yes
	 Equation = '{TKE_avg} = {TKE_avg} / |NumZones|'
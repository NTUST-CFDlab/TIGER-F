#!MC 1410
$!AlterData 
  IgnoreDivideByZero = Yes
  Equation = '{p_timeavg}= 0'
$! Loop |NumZones|
$! AlterData
   IgnoreDivideByZero = Yes
   Equation = '{p_timeavg} = {p_timeavg} +  {RHO}[|Loop|]'

$! EndLoop 
$! AlterData
     IgnoreDivideByZero = Yes
	 Equation = '{p_timeavg} = {p_timeavg} / |NumZones|'

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
	 
$!AlterData 
  IgnoreDivideByZero = Yes
  Equation = '{v_timeavg}= 0'
$! Loop |NumZones|
$! AlterData
   IgnoreDivideByZero = Yes
   Equation = '{v_timeavg} = {v_timeavg} +  {RHO-V}[|Loop|]'

$! EndLoop 
$! AlterData
     IgnoreDivideByZero = Yes
	 Equation = '{v_timeavg} = {v_timeavg} / |NumZones|'
	 
$!AlterData 
  IgnoreDivideByZero = Yes
  Equation = '{u_timeavg}= 0'
$! Loop |NumZones|
$! AlterData
   IgnoreDivideByZero = Yes
   Equation = '{u_timeavg} = {u_timeavg} +  {RHO-U}[|Loop|]'

$! EndLoop 
$! AlterData
     IgnoreDivideByZero = Yes
	 Equation = '{u_timeavg} = {u_timeavg} / |NumZones|'

$!AlterData 
  IgnoreDivideByZero = Yes
  Equation = '{kin_timeavg}= 0'
$! Loop |NumZones|
$! AlterData
   IgnoreDivideByZero = Yes
   Equation = '{kin_timeavg} = {kin_timeavg} +  {kin}[|Loop|]'

$! EndLoop 
$! AlterData
     IgnoreDivideByZero = Yes
	 Equation = '{kin_timeavg} = {kin_timeavg} / |NumZones|'
	 
$!AlterData 
  IgnoreDivideByZero = Yes
  Equation = '{wv_rss}= ({RHO-W}-{w_timeavg})*({RHO-V}-{v_timeavg})'
  
$!AlterData 
  IgnoreDivideByZero = Yes
  Equation = '{wv_rss_timeavg}= 0'
$! Loop |NumZones|
$! AlterData
   IgnoreDivideByZero = Yes
   Equation = '{wv_rss_timeavg} = {wv_rss_timeavg} +  {wv_rss}[|Loop|]'

$! EndLoop
$! AlterData
     IgnoreDivideByZero = Yes
	 Equation = '{wv_rss_timeavg} = {wv_rss_timeavg} / |NumZones|'
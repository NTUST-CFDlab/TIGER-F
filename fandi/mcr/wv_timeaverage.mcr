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
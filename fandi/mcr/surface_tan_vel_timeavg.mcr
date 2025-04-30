#!MC 1410
$!AlterData 
  IgnoreDivideByZero = Yes
  Equation = '{Tan_vel_avg}= 0'
$! Loop |NumZones|
$! AlterData
   IgnoreDivideByZero = Yes
   Equation = '{Tan_vel_avg} = {Tan_vel_avg} +  {Tan_vel}[|Loop|]'

$! EndLoop 
$! AlterData
     IgnoreDivideByZero = Yes
	 Equation = '{Tan_vel_avg} = {Tan_vel_avg} / |NumZones|'
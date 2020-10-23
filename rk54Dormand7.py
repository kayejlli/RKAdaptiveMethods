#  A family of embedded Runge-Kutta formulae, by J. R. Dormand and P. J. Prince, 
# Journal of Computational and Applied Mathematics, Vol 6, No. 1, 1980, pages 19 to 26.
# downloaded from http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK5/RKcoeff5d_1.pdf



k0=0.0E0
k1=0.2E0
k2=0.3E0
k3=0.8E0
k4=0.8888888888888888888888888888888888888888888888888888888888888888888888888888888888889E0
k5=1.0E0
k6=1.0E0
a1_0=0.2E0
a2_0=7.5E-2
a2_1=0.225E0
a3_0=0.9777777777777777777777777777777777777777777777777777777777777777777777777777777777778E0
a3_1=-3.733333333333333333333333333333333333333333333333333333333333333333333333333333333333E0
a3_2=3.555555555555555555555555555555555555555555555555555555555555555555555555555555555556E0
a4_0=2.952598689224203627495808565767413504039018442310623380582228318853833257125438195397E0
a4_1=-11.59579332418838591678097850937357110196616369455875628715134887974394147233653406493E0
a4_2=9.822892851699436061575979271452522481329065691205608901082152110958695320835238530712E0
a4_3=-0.2908093278463648834019204389574759945130315500685871056241426611796982167352537722908E0
a5_0=2.846275252525252525252525252525252525252525252525252525252525252525252525252525252525E0
a5_1=-10.75757575757575757575757575757575757575757575757575757575757575757575757575757575758E0
a5_2=8.906422717743472460453592529064227177434724604535925290642271774347246045359252906423E0
a5_3=0.2784090909090909090909090909090909090909090909090909090909090909090909090909090909091E0
a5_4=-0.2735313036020583190394511149228130360205831903945111492281303602058319039451114922813E0
a6_0=9.114583333333333333333333333333333333333333333333333333333333333333333333333333333333E-2
a6_1=0.0E0
a6_2=0.4492362982929020664869721473495058400718778077268643306379155435759209344115004492363E0
a6_3=0.6510416666666666666666666666666666666666666666666666666666666666666666666666666666667E0
a6_4=-0.3223761792452830188679245283018867924528301886792452830188679245283018867924528301887E0
a6_5=0.1309523809523809523809523809523809523809523809523809523809523809523809523809523809524E0
b0=9.114583333333333333333333333333333333333333333333333333333333333333333333333333333333E-2
b1=0.0E0
b2=0.4492362982929020664869721473495058400718778077268643306379155435759209344115004492363E0
b3=0.6510416666666666666666666666666666666666666666666666666666666666666666666666666666667E0
b4=-0.3223761792452830188679245283018867924528301886792452830188679245283018867924528301887E0
b5=0.1309523809523809523809523809523809523809523809523809523809523809523809523809523809524E0
b6=0.0E0
b_0=8.991319444444444444444444444444444444444444444444444444444444444444444444444444444444E-2
b_1=0.0E0
b_2=0.4534890685834082060497154836777478286912249176400119796346211440551063192572626534891E0
b_3=0.6140625E0
b_4=-0.2715123820754716981132075471698113207547169811320754716981132075471698113207547169811E0
b_5=8.904761904761904761904761904761904761904761904761904761904761904761904761904761904762E-2
b_6=2.5E-2



from CheckCoef import *
Print2(7, globals())
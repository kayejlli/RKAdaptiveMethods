MODULE rk1210FeaginMod

USE GlobalCommonMod
USE DyDtMod

IMPLICIT NONE
PRIVATE

!  Define access to SUBROUTINEs.

PUBLIC :: rk1210FeaginEachStep
! THE COEFFICIENTS OF RK12(10) TO 60 DIGITS
! using the notation of Fehlberg, Bettis, Horn, et alia


! k[i] are the nodes 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 k0=0.0D0,&
 k1=0.2D0,&
 k2=0.555555555555555555555555555555555555555555555555555555555556D0,&
 k3=0.833333333333333333333333333333333333333333333333333333333333D0,&
 k4=0.333333333333333333333333333333333333333333333333333333333333D0,&
 k5=1.0D0,&
 k6=0.671835709170513812712245661002797570438953420568682550710222D0,&
 k7=0.288724941110620201935458488967024976908118598341806976469674D0,&
 k8=0.5625D0,&
 k9=0.833333333333333333333333333333333333333333333333333333333333D0,&
 k10=0.947695431179199287562380162101836721649589325892740646458322D0,&
 k11=5.48112876863802643887753674810754475842153612931128785028369D-2,&
 k12=8.48880518607165350639838930162674302064148175640019542045934D-2,&
 k13=0.265575603264642893098114059045616835297201264164077621448665D0,&
 k14=0.5D0,&
 k15=0.734424396735357106901885940954383164702798735835922378551335D0,&
 k16=0.915111948139283464936016106983732569793585182435998045795407D0,&
 k17=0.947695431179199287562380162101836721649589325892740646458322D0,&
 k18=0.833333333333333333333333333333333333333333333333333333333333D0,&
 k19=0.288724941110620201935458488967024976908118598341806976469674D0,&
 k20=0.671835709170513812712245661002797570438953420568682550710222D0,&
 k21=0.333333333333333333333333333333333333333333333333333333333333D0,&
 k22=0.555555555555555555555555555555555555555555555555555555555556D0,&
 k23=0.2D0,&
 k24=1.0D0
! b[i] are coefficients for nodes k[i] 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 b0=2.38095238095238095238095238095238095238095238095238095238095D-2,&
 b1=2.34375D-2,&
 b2=3.125D-2,&
 b3=0.0D0,&
 b4=4.16666666666666666666666666666666666666666666666666666666667D-2,&
 b5=0.0D0,&
 b6=5.0D-2,&
 b7=5.0D-2,&
 b8=0.0D0,&
 b9=0.1D0,&
 b10=7.14285714285714285714285714285714285714285714285714285714286D-2,&
 b11=0.0D0,&
 b12=0.138413023680782974005350203145033146748813640089941234591267D0,&
 b13=0.215872690604931311708935511140681138965472074195773051123019D0,&
 b14=0.24380952380952380952380952380952380952380952380952380952381D0,&
 b15=0.215872690604931311708935511140681138965472074195773051123019D0,&
 b16=0.138413023680782974005350203145033146748813640089941234591267D0,&
 b17=-7.14285714285714285714285714285714285714285714285714285714286D-2,&
 b18=-0.1D0,&
 b19=-5.0D-2,&
 b20=-5.0D-2,&
 b21=-4.16666666666666666666666666666666666666666666666666666666667D-2,&
 b22=-3.125D-2,&
 b23=-2.34375D-2,&
 b24=2.38095238095238095238095238095238095238095238095238095238095D-2
! a[i,j] 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 a1_0=0.2D0,&
 a2_0=-0.216049382716049382716049382716049382716049382716049382716049D0,&
 a2_1=0.771604938271604938271604938271604938271604938271604938271605D0,&
 a3_0=0.208333333333333333333333333333333333333333333333333333333333D0,&
 a3_1=0.0D0,&
 a3_2=0.625D0,&
 a4_0=0.193333333333333333333333333333333333333333333333333333333333D0,&
 a4_1=0.0D0,&
 a4_2=0.22D0,&
 a4_3=-8.0D-2,&
 a5_0=0.1D0,&
 a5_1=0.0D0,&
 a5_2=0.0D0,&
 a5_3=0.4D0,&
 a5_4=0.5D0,&
 a6_0=0.103364471650010477570395435690481791543342708330349879244197D0,&
 a6_1=0.0D0,&
 a6_2=0.0D0,&
 a6_3=0.124053094528946761061581889237115328211074784955180298044074D0,&
 a6_4=0.48317116756103289928883648045196250872410925751728917730238D0,&
 a6_5=-3.8753024569476325208568144376762058039573330234136803880429D-2,&
 a7_0=0.12403826143183332408190458598017516814002467069863361229248D0,&
 a7_1=0.0D0,&
 a7_2=0.0D0,&
 a7_3=0.0D0,&
 a7_4=0.217050632197958486317846256953159942875916353757734167684657D0,&
 a7_5=1.37455792075966759812907801835048190594443990939408530842918D-2,&
 a7_6=-6.61095317267682844455831341498149531672668252085016565917546D-2,&
 a8_0=9.14774894856882983144991846980432197088832099976660100090486D-2,&
 a8_1=0.0D0,&
 a8_2=0.0D0,&
 a8_3=0.0D0,&
 a8_4=0.0D0,&
 a8_5=-5.44348523717469689965754944144838611346156873847009178068318D-3,&
 a8_6=6.80716801688453518578515120895103863112751730758794372203952D-2,&
 a8_7=0.408394315582641046727306852653894780093303185664924644551239D0,&
 a9_0=8.90013652502551018954509355423841780143232697403434118692699D-2,&
 a9_1=0.0D0,&
 a9_2=0.0D0,&
 a9_3=0.0D0,&
 a9_4=0.0D0,&
 a9_5=4.9952822664553236019779340842069280040589114940681409195581D-3,&
 a9_6=0.397918238819828997341739603001347156083435060931424970826304D0,&
 a9_7=0.427930210752576611068192608300897981558240730580396406312359D0,&
 a9_8=-8.65117637557827005740277475955029103267246394128995965941585D-2,&
 a10_0=6.95087624134907543112693906409809822706021061685544615255758D-2,&
 a10_1=0.0D0,&
 a10_2=0.0D0,&
 a10_3=0.0D0,&
 a10_4=0.0D0,&
 a10_5=0.129146941900176461970759579482746551122871751501482634045487D0,&
 a10_6=1.53073638102311295076342566143214939031177504112433874313011D0,&
 a10_7=0.577874761129140052546751349454576715334892100418571882718036D0,&
 a10_8=-0.951294772321088980532340837388859453930924498799228648050949D0,&
 a10_9=-0.408276642965631951497484981519757463459627174520978426909934D0,&
 a11_0=4.44861403295135866269453507092463581620165501018684152933313D-2,&
 a11_1=0.0D0,&
 a11_2=0.0D0,&
 a11_3=0.0D0,&
 a11_4=0.0D0,&
 a11_5=-3.80476867056961731984232686574547203016331563626856065717964D-3,&
 a11_6=1.06955064029624200721262602809059154469206077644957399593972D-2,&
 a11_7=2.09616244499904333296674205928919920806734650660039898074652D-2,&
 a11_8=-2.33146023259321786648561431551978077665337818756053603898847D-2,&
 a11_9=2.63265981064536974369934736325334761174975280887405725010964D-3,&
 a11_10=3.15472768977025060103545855572111407955208306374459723959783D-3,&
 a12_0=1.94588815119755475588801096525317761242073762016273186231215D-2,&
 a12_1=0.0D0,&
 a12_2=0.0D0,&
 a12_3=0.0D0,&
 a12_4=0.0D0,&
 a12_5=0.0D0,&
 a12_6=0.0D0,&
 a12_7=0.0D0,&
 a12_8=6.78512949171812509306121653452367476194364781259165332321534D-5,&
 a12_9=-4.2979585904927362327100533023016234356886338772488360367555D-5,&
 a12_10=1.76358982260285155407485928953302139937553442829975734148981D-5,&
 a12_11=6.53866627415027051009595231385181033549511358787382098351924D-2,&
 a13_0=0.206836835664277105916828174798272361078909196043446411598231D0,&
 a13_1=0.0D0,&
 a13_2=0.0D0,&
 a13_3=0.0D0,&
 a13_4=0.0D0,&
 a13_5=0.0D0,&
 a13_6=0.0D0,&
 a13_7=0.0D0,&
 a13_8=1.66796067104156472828045866664696450306326505094792505215514D-2,&
 a13_9=-8.79501563200710214457024178249986591130234990219959208704979D-3,&
 a13_10=3.46675455362463910824462315246379209427513654098596403637231D-3,&
 a13_11=-0.861264460105717678161432562258351242030270498966891201799225D0,&
 a13_12=0.908651882074050281096239478469262145034957129939256789178785D0,&
 a14_0=2.03926084654484010091511314676925686038504449562413004562382D-2,&
 a14_1=0.0D0,&
 a14_2=0.0D0,&
 a14_3=0.0D0,&
 a14_4=0.0D0,&
 a14_5=0.0D0,&
 a14_6=0.0D0,&
 a14_7=0.0D0,&
 a14_8=8.69469392016685948675400555583947505833954460930940959577347D-2,&
 a14_9=-1.91649630410149842286436611791405053287170076602337673587681D-2,&
 a14_10=6.55629159493663287364871573244244516034828755253746024098838D-3,&
 a14_11=9.87476128127434780903798528674033899738924968006632201445462D-2,&
 a14_12=5.35364695524996055083260173615567408717110247274021056118319D-3,&
 a14_13=0.301167864010967916837091303817051676920059229784957479998077D0,&
 a15_0=0.228410433917778099547115412893004398779136994596948545722283D0,&
 a15_1=0.0D0,&
 a15_2=0.0D0,&
 a15_3=0.0D0,&
 a15_4=0.0D0,&
 a15_5=0.0D0,&
 a15_6=0.0D0,&
 a15_7=0.0D0,&
 a15_8=-0.498707400793025250635016567442511512138603770959682292383042D0,&
 a15_9=0.134841168335724478552596703792570104791700727205981058201689D0,&
 a15_10=-3.8745824405583415843990422692402923093516105914280680567436D-2,&
 a15_11=-1.27473257473474844240388430824908952380979292713250350199641D0,&
 a15_12=1.43916364462877165201184452437038081875299303577911839630524D0,&
 a15_13=-0.214007467967990254219503540827349569639028092344812795499026D0,&
 a15_14=0.958202417754430239892724139109781371059908874605153648768037D0,&
 a16_0=2.00222477655974203614249646012506747121440306225711721209798D0,&
 a16_1=0.0D0,&
 a16_2=0.0D0,&
 a16_3=0.0D0,&
 a16_4=0.0D0,&
 a16_5=0.0D0,&
 a16_6=0.0D0,&
 a16_7=0.0D0,&
 a16_8=2.06701809961524912091954656438138595825411859673341600679555D0,&
 a16_9=0.62397813608613954195747127983149446615529231616702108066314D0,&
 a16_10=-4.62283685500311430283203554129062069391947101880112723185773D-2,&
 a16_11=-8.8497328836264961486007524672711894928660483545709270109463D0,&
 a16_12=7.74257707850855976227437225791835589560188590785037197433615D0,&
 a16_13=-0.588358519250869210993353314127711745644125882130941202896436D0,&
 a16_14=-1.10683733362380649395704708016953056176195769617014899442903D0,&
 a16_15=-0.929529037579203999778397238291233214220788057511899747507074D0,&
 a17_0=3.1378953341207344293445160898988879680816125933032210026831D0,&
 a17_1=0.0D0,&
 a17_2=0.0D0,&
 a17_3=0.0D0,&
 a17_4=0.0D0,&
 a17_5=0.129146941900176461970759579482746551122871751501482634045487D0,&
 a17_6=1.53073638102311295076342566143214939031177504112433874313011D0,&
 a17_7=0.577874761129140052546751349454576715334892100418571882718036D0,&
 a17_8=5.42088263055126683050056840891857421941300558851862156403363D0,&
 a17_9=0.231546926034829304872663800877643660904880180835945693836936D0,&
 a17_10=7.59292995578913560162301311785251873561801342333194895292058D-2,&
 a17_11=-12.3729973380186513287414553402595806591349822617535905976253D0,&
 a17_12=9.85455883464769543935957209317369202080367765721777101906955D0,&
 a17_13=8.5911143137043652957935770905236777288998049512232960115954D-2,&
 a17_14=-5.65242752862643921117182090081762761180392602644189218673969D0,&
 a17_15=-1.94300935242819610883833776782364287728724899124166920477873D0,&
 a17_16=-0.128352601849404542018428714319344620742146491335612353559923D0,&
 a18_0=1.38360054432196014878538118298167716825163268489922519995564D0,&
 a18_1=0.0D0,&
 a18_2=0.0D0,&
 a18_3=0.0D0,&
 a18_4=0.0D0,&
 a18_5=4.9952822664553236019779340842069280040589114940681409195581D-3,&
 a18_6=0.397918238819828997341739603001347156083435060931424970826304D0,&
 a18_7=0.427930210752576611068192608300897981558240730580396406312359D0,&
 a18_8=-1.30299107424475770916551439123047573342071475998399645982146D0,&
 a18_9=0.661292278669377029097112528107513072734573412294008071500699D0,&
 a18_10=-0.144559774306954349765969393688703463900585822441545655530145D0,&
 a18_11=-6.96576034731798203467853867461083919356792248105919255460819D0,&
 a18_12=6.65808543235991748353408295542210450632193197576935120716437D0,&
 a18_13=-1.66997375108841486404695805725510845049807969199236227575796D0,&
 a18_14=2.06413702318035263832289040301832647130604651223986452170089D0,&
 a18_15=-0.674743962644306471862958129570837723192079875998405058648892D0,&
 a18_16=-1.15618834794939500490703608435907610059605754935305582045729D-3,&
 a18_17=-5.44057908677007389319819914241631024660726585015012485938593D-3,&
 a19_0=0.951236297048287669474637975894973552166903378983475425758226D0,&
 a19_1=0.0D0,&
 a19_2=0.0D0,&
 a19_3=0.0D0,&
 a19_4=0.217050632197958486317846256953159942875916353757734167684657D0,&
 a19_5=1.37455792075966759812907801835048190594443990939408530842918D-2,&
 a19_6=-6.61095317267682844455831341498149531672668252085016565917546D-2,&
 a19_7=0.0D0,&
 a19_8=0.152281696736414447136604697040747131921486432699422112099617D0,&
 a19_9=-0.33774101835759984080230079313399800435464342445753966767008D0,&
 a19_10=-1.92825981633995781534949199286824400469353110630787982121133D-2,&
 a19_11=-3.68259269696866809932409015535499603576312120746888880201882D0,&
 a19_12=3.16197870406982063541533528419683854018352080342887002331312D0,&
 a19_13=-0.370462522106885290716991856022051125477943482284080569177386D0,&
 a19_14=-5.14974200365440434996434456698127984941168616474316871020314D-2,&
 a19_15=-8.29625532120152946787043541792848416659382675202720677536554D-4,&
 a19_16=2.79801041419278598986586589070027583961355402640879503213503D-6,&
 a19_17=4.18603916412360287969841020776788461794119440689356178942252D-2,&
 a19_18=0.279084255090877355915660874555379649966282167560126269290222D0,&
 a20_0=0.103364471650010477570395435690481791543342708330349879244197D0,&
 a20_1=0.0D0,&
 a20_2=0.0D0,&
 a20_3=0.124053094528946761061581889237115328211074784955180298044074D0,&
 a20_4=0.48317116756103289928883648045196250872410925751728917730238D0,&
 a20_5=-3.8753024569476325208568144376762058039573330234136803880429D-2,&
 a20_6=0.0D0,&
 a20_7=-0.438313820361122420391059788940960176420682836652600698580091D0,&
 a20_8=0.0D0,&
 a20_9=-0.218636633721676647685111485017151199362509373698288330593486D0,&
 a20_10=-3.12334764394719229981634995206440349766174759626578122323015D-2,&
 a20_11=0.0D0,&
 a20_12=0.0D0,&
 a20_13=0.0D0,&
 a20_14=0.0D0,&
 a20_15=0.0D0,&
 a20_16=0.0D0,&
 a20_17=3.12334764394719229981634995206440349766174759626578122323015D-2,&
 a20_18=0.218636633721676647685111485017151199362509373698288330593486D0,&
 a20_19=0.438313820361122420391059788940960176420682836652600698580091D0,&
 a21_0=0.193333333333333333333333333333333333333333333333333333333333D0,&
 a21_1=0.0D0,&
 a21_2=0.22D0,&
 a21_3=-8.0D-2,&
 a21_4=0.0D0,&
 a21_5=0.0D0,&
 a21_6=9.84256130499315928152900286856048243348202521491288575952143D-2,&
 a21_7=-0.196410889223054653446526504390100417677539095340135532418849D0,&
 a21_8=0.0D0,&
 a21_9=0.436457930493068729391826122587949137609670676712525034763317D0,&
 a21_10=6.5261372167572109856037093980555569835054381070841471673027D-2,&
 a21_11=0.0D0,&
 a21_12=0.0D0,&
 a21_13=0.0D0,&
 a21_14=0.0D0,&
 a21_15=0.0D0,&
 a21_16=0.0D0,&
 a21_17=-6.5261372167572109856037093980555569835054381070841471673027D-2,&
 a21_18=-0.436457930493068729391826122587949137609670676712525034763317D0,&
 a21_19=0.196410889223054653446526504390100417677539095340135532418849D0,&
 a21_20=-9.84256130499315928152900286856048243348202521491288575952143D-2,&
 a22_0=-0.216049382716049382716049382716049382716049382716049382716049D0,&
 a22_1=0.771604938271604938271604938271604938271604938271604938271605D0,&
 a22_2=0.0D0,&
 a22_3=0.0D0,&
 a22_4=-0.666666666666666666666666666666666666666666666666666666666667D0,&
 a22_5=0.0D0,&
 a22_6=-0.390696469295978451446999802258495981249099665294395945559163D0,&
 a22_7=0.0D0,&
 a22_8=0.0D0,&
 a22_9=0.0D0,&
 a22_10=0.0D0,&
 a22_11=0.0D0,&
 a22_12=0.0D0,&
 a22_13=0.0D0,&
 a22_14=0.0D0,&
 a22_15=0.0D0,&
 a22_16=0.0D0,&
 a22_17=0.0D0,&
 a22_18=0.0D0,&
 a22_19=0.0D0,&
 a22_20=0.390696469295978451446999802258495981249099665294395945559163D0,&
 a22_21=0.666666666666666666666666666666666666666666666666666666666667D0,&
 a23_0=0.2D0,&
 a23_1=0.0D0,&
 a23_2=-0.164609053497942386831275720164609053497942386831275720164609D0,&
 a23_3=0.0D0,&
 a23_4=0.0D0,&
 a23_5=0.0D0,&
 a23_6=0.0D0,&
 a23_7=0.0D0,&
 a23_8=0.0D0,&
 a23_9=0.0D0,&
 a23_10=0.0D0,&
 a23_11=0.0D0,&
 a23_12=0.0D0,&
 a23_13=0.0D0,&
 a23_14=0.0D0,&
 a23_15=0.0D0,&
 a23_16=0.0D0,&
 a23_17=0.0D0,&
 a23_18=0.0D0,&
 a23_19=0.0D0,&
 a23_20=0.0D0,&
 a23_21=0.0D0,&
 a23_22=0.164609053497942386831275720164609053497942386831275720164609D0,&
 a24_0=1.47178724881110408452949550989023611293535315518571691939396D0,&
 a24_1=0.7875D0,&
 a24_2=0.421296296296296296296296296296296296296296296296296296296296D0,&
 a24_3=0.0D0,&
 a24_4=0.291666666666666666666666666666666666666666666666666666666667D0,&
 a24_5=0.0D0,&
 a24_6=0.348600717628329563206854421629657569274689947367847465753757D0,&
 a24_7=0.22949954476899484958289023371055544707382356966650670066251D0,&
 a24_8=5.79046485790481979159831978177003471098279506036722411333192D0,&
 a24_9=0.418587511856506868874073759426596207226461447604248151080016D0,&
 a24_10=0.307039880222474002649653817490106690389251482313213999386651D0,&
 a24_11=-4.68700905350603332214256344683853248065574415794742040470287D0,&
 a24_12=3.13571665593802262152038152399873856554395436199962915429076D0,&
 a24_13=1.40134829710965720817510506275620441055845017313930508348898D0,&
 a24_14=-5.5293110143949902362901030600576433642127605577765815640091D0,&
 a24_15=-0.853138235508063349309546894974784906188927508039552519557498D0,&
 a24_16=0.103575780373610140411804607167772795518293914458500175573749D0,&
 a24_17=-0.140474416950600941142546901202132534870665923700034957196546D0,&
 a24_18=-0.418587511856506868874073759426596207226461447604248151080016D0,&
 a24_19=-0.22949954476899484958289023371055544707382356966650670066251D0,&
 a24_20=-0.348600717628329563206854421629657569274689947367847465753757D0,&
 a24_21=-0.291666666666666666666666666666666666666666666666666666666667D0,&
 a24_22=-0.421296296296296296296296296296296296296296296296296296296296D0,&
 a24_23=-0.7875D0
! The estimate of the local truncation error is  (49/640)  h ( f(t1,x1)-f(t23,x23) )


CONTAINS

SUBROUTINE rk1210FeaginEachStep(y0,yn,h,hnew,rerun,test)
 REAL(KIND=8), DIMENSION(:), INTENT(IN) :: y0     ! y(t)
 REAL(KIND=8), DIMENSION(SIZE(y0)), INTENT(OUT) :: yn ! y(t+h)
 REAL(KIND=8), INTENT(IN) :: h ! initial step size
 REAL(KIND=8), INTENT(OUT) :: hnew       ! new step size
 LOGICAL, INTENT(OUT) :: test, rerun ! test=stop the program, rerun=re-run this step, reject
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: tolh ! the tolerance, determined by the problem
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: yerr ! the error between embedded method
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: ynp  ! the embedded
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: ymax ! the max value among y0 and yn
 REAL(KIND=8) :: err
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: y1,y2,y3,y4,y5,y6,y7,y8,y9, &
                                      y10,y11,y12,y13,y14,y15,y16,y17,y18, &
                                      y19,y20,y21,y22,y23,y24
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: dy0,dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8, &
                                      dy9,dy10,dy11,dy12,dy13,dy14,dy15,dy16,dy17, &
                                      dy18,dy19,dy20,dy21,dy22,dy23,dy24
 INTEGER :: i
  test = .False.
  rerun = .False. 
  ! use y0 to get dy0
  CALL  dev(y0,dy0,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y1=y0+h*(a1_0*dy0)
  ! use y1 to get dy1
  CALL  dev(y1,dy1,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y2=y0+h*(a2_0*dy0+a2_1*dy1)
  ! use y2 to get dy2
  CALL  dev(y2,dy2,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y3=y0+h*(a3_0*dy0+a3_1*dy1+a3_2*dy2)
  ! use y3 to get dy3
  CALL  dev(y3,dy3,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y4=y0+h*(a4_0*dy0+a4_1*dy1+a4_2*dy2+a4_3*dy3)
  ! use y4 to get dy4
  CALL  dev(y4,dy4,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y5=y0+h*(a5_0*dy0+a5_1*dy1+a5_2*dy2+a5_3*dy3+a5_4*dy4)
  ! use y5 to get dy5
  CALL  dev(y5,dy5,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y6=y0+h*(a6_0*dy0+a6_1*dy1+a6_2*dy2+a6_3*dy3+a6_4*dy4+a6_5*dy5)
  ! use y6 to get dy6
  CALL  dev(y6,dy6,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y7=y0+h*(a7_0*dy0+a7_1*dy1+a7_2*dy2+a7_3*dy3+a7_4*dy4+a7_5*dy5 + &
        &  a7_6*dy6)
  ! use y7 to get dy7
  CALL  dev(y7,dy7,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y8=y0+h*(a8_0*dy0+a8_1*dy1+a8_2*dy2+a8_3*dy3+a8_4*dy4+a8_5*dy5 + &
        &  a8_6*dy6+a8_7*dy7)
  ! use y8 to get dy8
  CALL  dev(y8,dy8,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y9=y0+h*(a9_0*dy0+a9_1*dy1+a9_2*dy2+a9_3*dy3+a9_4*dy4+a9_5*dy5 + &
        &  a9_6*dy6+a9_7*dy7+a9_8*dy8)
  ! use y9 to get dy9
  CALL  dev(y9,dy9,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y10=y0+h*(a10_0*dy0+a10_1*dy1+a10_2*dy2+a10_3*dy3+a10_4*dy4+a10_5*dy5 + &
         &  a10_6*dy6+a10_7*dy7+a10_8*dy8+a10_9*dy9)
  ! use y10 to get dy10
  CALL  dev(y10,dy10,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y11=y0+h*(a11_0*dy0+a11_1*dy1+a11_2*dy2+a11_3*dy3+a11_4*dy4+a11_5*dy5 + &
         &  a11_6*dy6+a11_7*dy7+a11_8*dy8+a11_9*dy9+a11_10*dy10)
  ! use y11 to get dy11
  CALL  dev(y11,dy11,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y12=y0+h*(a12_0*dy0+a12_1*dy1+a12_2*dy2+a12_3*dy3+a12_4*dy4+a12_5*dy5 + &
         &  a12_6*dy6+a12_7*dy7+a12_8*dy8+a12_9*dy9+a12_10*dy10+a12_11*dy11)
  ! use y12 to get dy12
  CALL  dev(y12,dy12,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y13=y0+h*(a13_0*dy0+a13_1*dy1+a13_2*dy2+a13_3*dy3+a13_4*dy4+a13_5*dy5 + &
         &  a13_6*dy6+a13_7*dy7+a13_8*dy8+a13_9*dy9+a13_10*dy10+a13_11*dy11 + &
         &  a13_12*dy12)
  ! use y13 to get dy13
  CALL  dev(y13,dy13,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y14=y0+h*(a14_0*dy0+a14_1*dy1+a14_2*dy2+a14_3*dy3+a14_4*dy4+a14_5*dy5 + &
         &  a14_6*dy6+a14_7*dy7+a14_8*dy8+a14_9*dy9+a14_10*dy10+a14_11*dy11 + &
         &  a14_12*dy12+a14_13*dy13)
  ! use y14 to get dy14
  CALL  dev(y14,dy14,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y15=y0+h*(a15_0*dy0+a15_1*dy1+a15_2*dy2+a15_3*dy3+a15_4*dy4+a15_5*dy5 + &
         &  a15_6*dy6+a15_7*dy7+a15_8*dy8+a15_9*dy9+a15_10*dy10+a15_11*dy11 + &
         &  a15_12*dy12+a15_13*dy13+a15_14*dy14)
  ! use y15 to get dy15
  CALL  dev(y15,dy15,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y16=y0+h*(a16_0*dy0+a16_1*dy1+a16_2*dy2+a16_3*dy3+a16_4*dy4+a16_5*dy5 + &
         &  a16_6*dy6+a16_7*dy7+a16_8*dy8+a16_9*dy9+a16_10*dy10+a16_11*dy11 + &
         &  a16_12*dy12+a16_13*dy13+a16_14*dy14+a16_15*dy15)
  ! use y16 to get dy16
  CALL  dev(y16,dy16,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y17=y0+h*(a17_0*dy0+a17_1*dy1+a17_2*dy2+a17_3*dy3+a17_4*dy4+a17_5*dy5 + &
         &  a17_6*dy6+a17_7*dy7+a17_8*dy8+a17_9*dy9+a17_10*dy10+a17_11*dy11 + &
         &  a17_12*dy12+a17_13*dy13+a17_14*dy14+a17_15*dy15+a17_16*dy16)
  ! use y17 to get dy17
  CALL  dev(y17,dy17,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y18=y0+h*(a18_0*dy0+a18_1*dy1+a18_2*dy2+a18_3*dy3+a18_4*dy4+a18_5*dy5 + &
         &  a18_6*dy6+a18_7*dy7+a18_8*dy8+a18_9*dy9+a18_10*dy10+a18_11*dy11 + &
         &  a18_12*dy12+a18_13*dy13+a18_14*dy14+a18_15*dy15+a18_16*dy16+a18_17*dy17)
  ! use y18 to get dy18
  CALL  dev(y18,dy18,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y19=y0+h*(a19_0*dy0+a19_1*dy1+a19_2*dy2+a19_3*dy3+a19_4*dy4+a19_5*dy5 + &
         &  a19_6*dy6+a19_7*dy7+a19_8*dy8+a19_9*dy9+a19_10*dy10+a19_11*dy11 + &
         &  a19_12*dy12+a19_13*dy13+a19_14*dy14+a19_15*dy15+a19_16*dy16+a19_17*dy17 + &
         &  a19_18*dy18)
  ! use y19 to get dy19
  CALL  dev(y19,dy19,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y20=y0+h*(a20_0*dy0+a20_1*dy1+a20_2*dy2+a20_3*dy3+a20_4*dy4+a20_5*dy5 + &
         &  a20_6*dy6+a20_7*dy7+a20_8*dy8+a20_9*dy9+a20_10*dy10+a20_11*dy11 + &
         &  a20_12*dy12+a20_13*dy13+a20_14*dy14+a20_15*dy15+a20_16*dy16+a20_17*dy17 + &
         &  a20_18*dy18+a20_19*dy19)
  ! use y20 to get dy20
  CALL  dev(y20,dy20,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y21=y0+h*(a21_0*dy0+a21_1*dy1+a21_2*dy2+a21_3*dy3+a21_4*dy4+a21_5*dy5 + &
         &  a21_6*dy6+a21_7*dy7+a21_8*dy8+a21_9*dy9+a21_10*dy10+a21_11*dy11 + &
         &  a21_12*dy12+a21_13*dy13+a21_14*dy14+a21_15*dy15+a21_16*dy16+a21_17*dy17 + &
         &  a21_18*dy18+a21_19*dy19+a21_20*dy20)
  ! use y21 to get dy21
  CALL  dev(y21,dy21,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y22=y0+h*(a22_0*dy0+a22_1*dy1+a22_2*dy2+a22_3*dy3+a22_4*dy4+a22_5*dy5 + &
         &  a22_6*dy6+a22_7*dy7+a22_8*dy8+a22_9*dy9+a22_10*dy10+a22_11*dy11 + &
         &  a22_12*dy12+a22_13*dy13+a22_14*dy14+a22_15*dy15+a22_16*dy16+a22_17*dy17 + &
         &  a22_18*dy18+a22_19*dy19+a22_20*dy20+a22_21*dy21)
  ! use y22 to get dy22
  CALL  dev(y22,dy22,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y23=y0+h*(a23_0*dy0+a23_1*dy1+a23_2*dy2+a23_3*dy3+a23_4*dy4+a23_5*dy5 + &
         &  a23_6*dy6+a23_7*dy7+a23_8*dy8+a23_9*dy9+a23_10*dy10+a23_11*dy11 + &
         &  a23_12*dy12+a23_13*dy13+a23_14*dy14+a23_15*dy15+a23_16*dy16+a23_17*dy17 + &
         &  a23_18*dy18+a23_19*dy19+a23_20*dy20+a23_21*dy21+a23_22*dy22)
  ! use y23 to get dy23
  CALL  dev(y23,dy23,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y24=y0+h*(a24_0*dy0+a24_1*dy1+a24_2*dy2+a24_3*dy3+a24_4*dy4+a24_5*dy5 + &
         &  a24_6*dy6+a24_7*dy7+a24_8*dy8+a24_9*dy9+a24_10*dy10+a24_11*dy11 + &
         &  a24_12*dy12+a24_13*dy13+a24_14*dy14+a24_15*dy15+a24_16*dy16+a24_17*dy17 + &
         &  a24_18*dy18+a24_19*dy19+a24_20*dy20+a24_21*dy21+a24_22*dy22+a24_23*dy23)
  ! use y24 to get dy24
  CALL  dev(y24,dy24,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  yn=y0+h*(b0*dy0+b1*dy1+b2*dy2+b3*dy3+b4*dy4+b5*dy5+b6*dy6+b7*dy7 + &
        &  b8*dy8+b9*dy9+b10*dy10+b11*dy11+b12*dy12+b13*dy13+b14*dy14+b15*dy15 + &
        &  b16*dy16+b17*dy17+b18*dy18+b19*dy19+b20*dy20+b21*dy21+b22*dy22+b23*dy23 + &
        &  b24*dy24)
  yerr = (4.9D1/6.4D2)*h*ABS(dy1-dy23)
  ! Find the max value of y among this step
  DO i = 1, SIZE(y0)
    ymax(i) = MAX(MAX(ABS(y0(i)), ABS(yn(i))), h*ABS(dy0(i)))
  END DO
  tolh = rtol*ymax + atol(1:SIZE(y0)) ! atol might be a longer array

  ! using the error to estimate the next step
  err = MAXVAL(ABS(yerr/tolh))
  IF (err.GT.1.D0) THEN
    rerun = .True.
    hnew = MAX(0.9D0*err**(-9.0909090909091D-02), 0.1D0)*h ! no less than factor of 0.1
    ! PRINT *, 'Decrease time step by', 0.9D0*err**(-9.0909090909091D-02),MAX(0.9D0*err**(-9.0909090909091D-02), 0.1D0)
  ELSE
    rerun = .False.
    hnew = MIN(5.D0, 0.9D0*err**(-9.0909090909091D-02))*h ! no more than factor of 5
    ! PRINT *, 'Increase time step by', 0.9D0*err**(-9.0909090909091D-02),MIN(5.D0,0.9D0*err**(-9.0909090909091D-02))
  END IF

  ! adjust the step
  hnew = MAX(MIN(hnew,MaxStepSize),MinStepSize)
  IF (rerun) THEN
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) THEN ! h is already the min
      test = .True. ! stop the program
    ELSE
      hnew = MAX(MinStepSize,MIN(hnew, h*0.5D0)) ! if have not reduced, reduce to half
    END IF
    RETURN
  END IF

  ! check if any value have went crazy (Nan or Inf)
  DO i = 1, SIZE(y0)
    CALL checkNanInf(yn(i), rerun)
    IF (rerun) THEN
      IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) THEN ! h is already the min
        test = .True. ! stop the program
      ELSE
        hnew = MAX(MinStepSize,MIN(hnew, h*0.5D0)) ! if have not reduced, reduce to half
      END IF
      RETURN
    END IF
  END DO

  RETURN
END SUBROUTINE

END MODULE
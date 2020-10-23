# http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK8/RKcoeff8d_1.pdf
import numpy as np

CopyRightInfo = '\
# Prince Dormand\n\
# High order embedded Runge-Kutta formulae\n\
# by P.J.Prince and J.R.Dormand, Journal of Computational and Applied Mathematics, vol. 7, 1981, pages 67-75.\n\
# downloaded from http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK8/RKcoeff8d_1.pdf'

print(CopyRightInfo) 
print('\n'*2)

a='\
c[1]=0|\
c[2]=1/18|\
c[3]=1/12|\
c[4]=1/8|\
c[5]=5/16|\
c[6]=3/8|\
c[7]=59/400|\
c[8]=93/200|\
c[9]=5490023248/9719169821|\
c[10]=13/20|\
c[11]=30992876149296355/33518267164510641|\
c[12]=1|\
c[13]=1|\
a[2,1]=1/18|\
a[3,1]=1/48|\
a[3,2]=1/16|\
a[4,1]=1/32|\
a[4,2]=0|\
a[4,3]=3/32|\
a[5,1]=5/16|\
a[5,2]=0|\
a[5,3]=-75/64|\
a[5,4]=75/64|\
a[6,1]=3/80|\
a[6,2]=0|\
a[6,3]=0|\
a[6,4]=3/16|\
a[6,5]=3/20|\
a[7,1]=215595617/4500000000|\
a[7,2]=0|\
a[7,3]=0|\
a[7,4]=202047683/1800000000|\
a[7,5]=-28693883/1125000000|\
a[7,6]=23124283/1800000000|\
a[8,1]=14873762658037143/879168438156250000|\
a[8,2]=0|\
a[8,3]=0|\
a[8,4]=3467633544794897/8940695981250000|\
a[8,5]=1474287494383247/40978189914062500|\
a[8,6]=26709270507070017/135600555715625000|\
a[8,7]=-14591655588284/84484570233063|\
a[9,1]=7586331039021946882049083502441337664277676907617750536566352/109794461601491217860220353338581031394059220336451160078730445|\
a[9,2]=0|\
a[9,3]=0|\
a[9,4]=-236057339412812449835946465344221735535939129430991059693568/372184615598275314780407977418918750488336340123563254504171|\
a[9,5]=-3299739166368883603096250588167927276977533790499480498577408/20470153857905142312922438758040531276858498706795978997729405|\
a[9,6]=4695919603694846215470554638065271273971468502369170235542016/33868800019443053645017125945121606294438606951244256159879561|\
a[9,7]=291851811898394201384602939640627532330843113837053004434432000000/310174233778061645620360730197195350622945922304711702829528117367|\
a[9,8]=6992959981041103840944260661352231159203510904000000/33042342481018810238716485165383193327572243242031481|\
a[10,1]=99299034813490800741867453179778547/540971123539151162906952826011200000|\
a[10,2]=0|\
a[10,3]=0|\
a[10,4]=-2493835259080554724582/1010153717930905426875|\
a[10,5]=-48550347897506146536052/166675363458599395434375|\
a[10,6]=-24871192635697392099560348960246/939492072180864357472739828818125|\
a[10,7]=478776089216929482237673925052922000/168119099731629344552415590032785027|\
a[10,8]=6560308981643238155096750/23314158982833116227901307|\
a[10,9]=1586281686644478270321241459439899956623408540189275177/12818966182821619734532382093543907143647820508227904000|\
a[11,1]=-102116003386322998978127600084904875522141269364791505043913504184525097818434721165778087547359160299919872547571820573487921693/84016717385376362440519288454722754561118206109968455863629915569413007015484884082989277327146750610032897620427741658059440288|\
a[11,2]=0|\
a[11,3]=0|\
a[11,4]=338590872606752219742507143357021902717271169524361004010718467428498066558752974165816979255870352236800/20308212073515087965058545521329962060416948491603802421256875704911573108931922671691153944392874968051|\
a[11,5]=68189290605616416787948548385820859588684790288743680764422049661712817412461535969920258787664375619072/74463444269555322538548000244876527554862144469213942211275210918009101399417049796200897796107208216187|\
a[11,6]=-1734282043732424474072631498514610096486835338935534079058145376760622893524503274680375038942168945756187943481380463560951840/286345537377499805912462279621622489249909215975695809863482134802066603511244489020404711919081949534640172065152437496912477|\
a[11,7]=-3399549280223124443696423490103003766707892326374755946138975000967466690241111348721006509128775254952212682658842765965521154240000000/212424385105117691648087703103838079790425456287912424851546922389328485500145214289225448961304538830766072442444722564103495915888123|\
a[11,8]=14452808190943733856347403293564049428070036006455540637351575894308889412108389906599600485253194980566957563315340127500000/973298753951638431793701721528200883789914680313298926814615071301495341142665245758696799918623095581715765886887649741383|\
a[11,9]=-847205714160239289113307424793539077951658318917591980262304042838612275700008766016957700930195545053374220841398660187944621107065829310608865394026418258355/63358704383980726998416112830322706485300332630289060627019459285960825979588560697460438306253611095891491565590971432387489415884103732012574255897878321824|\
a[11,10]=115188988949323598098458035263894669359112068207548636038131244599058496172710646622536373145562218909633738697549245770000/22435701423704647109276644681016984863989966659062291511947760084943925084166270812354794844590216383205333034660617626349|\
a[12,1]=21969012306961489525323859125985377266525845354279828748/84868015648089839210997460517819380601933600521692915045|\
a[12,2]=0|\
a[12,3]=0|\
a[12,4]=-2291872762438069505504/480025046760766258851|\
a[12,5]=-3829018311866050387904/8800459190614048078935|\
a[12,6]=-607977714773374460437401016185253441418120832060126402968/199370728929424959394190105343852509479613745231838418799|\
a[12,7]=5302029233035772894614097632213626682295966947853615180783170000000/950538766256052885387161080614691196420735587733978871061913292363|\
a[12,8]=102968047255116137164987219663037502898143843145000000/16726911019578511096352500731821705820659977305290973|\
a[12,9]=-111383789341965407321602142444917514115800834690201329379027449761759895100011973929185171163615/22003454775272439861723739055800175619777853128055268766511800511549546753240522083740083243539|\
a[12,10]=44737471541467333111555512048686345065750/20391511842264262870398286145868618178341|\
a[12,11]=596546910748352988538198147432444829112451075399436970876618894337461087953328002664759407401623072330633057948252/4431076125983762085449284205348478790535717302043416234911901479328512794465980800998816354448181196721636373483787|\
a[13,1]=1066221205855832326088695778460159015192405644968016897066521076847764032613686056268693633/1296431693610525557488309197474904206216262654240544950471874305723890174339356551609704000|\
a[13,2]=0|\
a[13,3]=0|\
a[13,4]=-1335791413506612664643690684478806471077526746614666064/114574907798601779179110271814903983120429559544320175|\
a[13,5]=-1591415543044168099882026495959288688569084060473110176/2100539976307699284950354983273239690541208591645869875|\
a[13,6]=33975758488532631832742416857645572913178866704247539610423012370193845167470455176890924/47586856225469573819304596274208152402640120925455970356063642741972959597009066064956075|\
a[13,7]=12176653428667113090492984656207574633063967759246601254930448409444470870786024235115138527800000/1008353786145118968620988891518234034224047994442049071310258686840184337101721351612973016221399|\
a[13,8]=-339784374935367314296824613776444883113869450234942131172912300100535979345925250000/159698690787587746004588725210359673189662237866695585709500421500486548151424426361|\
a[13,9]=4955095692700499418628052380948016677978733013841365878109775677669056866398110949788869771135857671298802131693154421086808143/2489789885462873158531234022579722982784822257458164105126884288597324542930882581099522281388970940826324647386340365850671680|\
a[13,10]=-563115171027780776675066866318087406247194110301648522108648094708415/2403532595444498372383116767918060257292523183751650851596520916634577|\
a[13,11]=147332487580158450887955957061658718012538967463083369806963200702426559434915876714751833908862217396388157664714990174448521780809/837599084085749358149340415048050308970085851893614803629073546048735327947816070400330404870816820234727495143522673498826476267825|\
a[13,12]=0|\
b[1]=212810988215683677989664967567559/5097575504458999984164528930580800|\
b[2]=0|\
b[3]=0|\
b[4]=0|\
b[5]=0|\
b[6]=-570667999368605802515460802224128/10291145812277763122885317774476825|\
b[7]=3970894643399159150754126826496000000000000/16592904867230933191457493387696939021741363|\
b[8]=177094288219480472437690862000000000000/251729356670100506734814442705774463449|\
b[9]=-66822609448295850920212176513645119787713273203022994500406050793972052314809461629969645683/87952305220338336969447643899150816363456821562985998778022435070001091778042097545895594560|\
b[10]=314652731163869955629145958568800000/476340207420551356675670184044905167|\
b[11]=177014954088789647707522848990757432519504314686067075784476503038212450536095365316360385634933688213244039743969578872631174179769/1119019983628991838522384101261104859676427163726922121733732080377576616485631933067985100908132443862205090961383250990215178108200|\
b[12]=-454665916000392064556420344242099/1909482158429176288068071462671400|\
b[13]=1/4|\
b_[1]=7136040226482108704342809557217/241464102842794736092004001974880|\
b_[2]=0|\
b_[3]=0|\
b_[4]=0|\
b_[5]=0|\
b_[6]=-15349154422148033115423212285265536/18524062462099973621193571994058285|\
b_[7]=45434521806506196832804182374790400000000/145978635195580057402851847985603569106229|\
b_[8]=365696286946774693155766999232150000000/148214481030059176862554298041717674741|\
b_[9]=-836336669851503831866889530158468123932231502753408325817124013619515886965077571/328368994730082689886153304749497093954319862912916225944630536728837081959128864|\
b_[10]=294694385044387823293019951454286000/204145803180236295718144364590673643|\
b_[11]=1759482754698187564675489259591170188433054767657805212470918093603353527288272972728828708146708084742711724049636/22155380629918810427246421026742393952678586510217081174559507396642563972329904004994081772240905983608181867418935|\
b_[12]=2/45|\
b_[13]=0|\
'

from PeterStone import * 
PrintFixed(a)
print('\n'*2)

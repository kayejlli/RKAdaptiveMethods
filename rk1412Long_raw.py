
CopyRightInfo = '\
# Terry Feagin 14/12 th order\n\
# downloaded from http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK14/RKcoeff14a_1.pdf' 


a='\
c[1]=0|\
c[2]=.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111|\
c[3]=.5555555555555555555555555555555555555555555555555555555555555555555555555555555555556|\
c[4]=.8333333333333333333333333333333333333333333333333333333333333333333333333333333333333|\
c[5]=.3333333333333333333333333333333333333333333333333333333333333333333333333333333333333|\
c[6]=1.|\
c[7]=.6699869792727729217646837855059985139388452296384603532851421391683474428303956826239|\
c[8]=.2970683842138183573895847168082194132233320946989156873791682903324708698499266217383|\
c[9]=.7272727272727272727272727272727272727272727272727272727272727272727272727272727272727|\
c[10]=.1401527990421887652761874879669467176298064630825329362873230163439023340348096838456|\
c[11]=.7007010397701507371510998548307493379414070492655464089692218490447945746638665522966|\
c[12]=.3636363636363636363636363636363636363636363636363636363636363636363636363636363636364|\
c[13]=.2631578947368421052631578947368421052631578947368421052631578947368421052631578947368|\
c[14]=.392172246650270859125196642501208648863714315266128052078483e-1|\
c[15]=.8129175029283767629833931592780365061896123726172385507744269795906758195776958783707|\
c[16]=.1666666666666666666666666666666666666666666666666666666666666666666666666666666666667|\
c[17]=.9|\
c[18]=.6412992574519669233127711938966828094810966516150832254029235721305050295351572963693e-1|\
c[19]=.2041499092834288489277446343010234050271495052413337516288702042649259099754335560687|\
c[20]=.3953503910487605656156713698273243723522272974566594505545766538389345381768585023057|\
c[21]=.6046496089512394343843286301726756276477727025433405494454233461610654618231414976943|\
c[22]=.7958500907165711510722553656989765949728504947586662483711297957350740900245664439313|\
c[23]=.9358700742548033076687228806103317190518903348384916774597076427869494970464842703631|\
c[24]=.1666666666666666666666666666666666666666666666666666666666666666666666666666666666667|\
c[25]=.8129175029283767629833931592780365061896123726172385507744269795906758195776958783707|\
c[26]=.392172246650270859125196642501208648863714315266128052078483e-1|\
c[27]=.3636363636363636363636363636363636363636363636363636363636363636363636363636363636364|\
c[28]=.7007010397701507371510998548307493379414070492655464089692218490447945746638665522966|\
c[29]=.1401527990421887652761874879669467176298064630825329362873230163439023340348096838456|\
c[30]=.2970683842138183573895847168082194132233320946989156873791682903324708698499266217383|\
c[31]=.6699869792727729217646837855059985139388452296384603532851421391683474428303956826239|\
c[32]=.3333333333333333333333333333333333333333333333333333333333333333333333333333333333333|\
c[33]=.5555555555555555555555555555555555555555555555555555555555555555555555555555555555556|\
c[34]=.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111|\
c[35]=1.|\
a[2,1]=.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111|\
a[3,1]=-.8333333333333333333333333333333333333333333333333333333333333333333333333333333333333|\
a[3,2]=1.388888888888888888888888888888888888888888888888888888888888888888888888888888888889|\
a[4,1]=.2083333333333333333333333333333333333333333333333333333333333333333333333333333333333|\
a[4,2]=0.|\
a[4,3]=.625|\
a[5,1]=.1933333333333333333333333333333333333333333333333333333333333333333333333333333333333|\
a[5,2]=0.|\
a[5,3]=.22|\
a[5,4]=-.8e-1|\
a[6,1]=.1|\
a[6,2]=0.|\
a[6,3]=0.|\
a[6,4]=.4|\
a[6,5]=.5|\
a[7,1]=.1034845616366797766729935465119103444997447982019713166066629728281981965079290745983|\
a[7,2]=0.|\
a[7,3]=0.|\
a[7,4]=.1220688873064072225896440828689620771395927148341621347412746563709055937325311521675|\
a[7,5]=.4825744903312466224751347801256881128659190238501680496794015023696413273862321544150|\
a[7,6]=-.3814096000156069997308862400056202056641130724784114774219699240039767479629669855696e-1|\
a[8,1]=.1243805266540944128815164208687993162684914663596714231632892354628068537117612942798|\
a[8,2]=0.|\
a[8,3]=0.|\
a[8,4]=0.|\
a[8,5]=.2261202821975843014222386629792029011967523207426331439651447460281196206643404356021|\
a[8,6]=.1378858876180808806076958370164778145309694174914933853635428709475288586061552782365e-1|\
a[8,7]=-.6722101339966844497493995074143058569500863415253821828561997825320849038679063596730e-1|\
a[9,1]=.9369190656596738155308854560830059338663496952177500856556033862893464429241815101000e-1|\
a[9,2]=0.|\
a[9,3]=0.|\
a[9,4]=0.|\
a[9,5]=0.|\
a[9,6]=-.6134068434505109872294989956416647356209145071288588710070986068372475355320835997035e-2|\
a[9,7]=.2160198256255030637088600976598665734909794332781173201886676706066128640340557614360|\
a[9,8]=.4236950635157619373376190739609767532058674695441235326831157041055522397561196508237|\
a[10,1]=.8384798124090526646169687913728140859805331392249111310693346670107922625197375034871e-1|\
a[10,2]=0.|\
a[10,3]=0.|\
a[10,4]=0.|\
a[10,5]=0.|\
a[10,6]=-.1179493671009738143197550560312957753679619605907361507776128268875265788248790903515e-1|\
a[10,7]=-.2472990205688126523394738387431945983259928403533401326974984247503501083158412965835|\
a[10,8]=.9780808583677290122593130140812916655037406554767339407565991037499621093437371932341e-1|\
a[10,9]=.2175906892434206313600086517678603183441681200247821768799893467069296630467914197921|\
a[11,1]=.6152553597694282279545623896143147143334239690648211074539397569215087099333654844097e-1|\
a[11,2]=0.|\
a[11,3]=0.|\
a[11,4]=0.|\
a[11,5]=0.|\
a[11,6]=.5922327803245033080429900057980465247383895604442571368349896773084347972825775455007e-2|\
a[11,7]=.4703261599638411122172243032058941134553625307461088250108483236601604516650193568134|\
a[11,8]=.2996888638486790008539818370961923991368311216717812791841936858888827504094204242461|\
a[11,9]=-.2476568775939949146899922763298108258539580692639470955481886317480090967647905771626|\
a[11,10]=.1108950297714376828939998518390617145224451736006787182086245987785252503880550245038|\
a[12,1]=.4197000733627825798617928647872777872134836565431046112459945389674655429048057710370e-1|\
a[12,2]=0.|\
a[12,3]=0.|\
a[12,4]=0.|\
a[12,5]=0.|\
a[12,6]=-.317987696266205093901912847692712407988609169703103952205634e-2|\
a[12,7]=.8063977149061920772608217115203795063935431115674197501197468839656405367779525213500|\
a[12,8]=.9759831264123889790935228506842888513146720480030545503571875185550549213299958241991e-1|\
a[12,9]=.7785755781583989090275124464529272389997634605941819649588520345133050850477185489203|\
a[12,10]=.2048904238315994281894992020981056033120292350814206535748293420400885242747823516625|\
a[12,11]=-1.562615796274681883070709439505278252114628922364243608928053762634922556160297217820|\
a[13,1]=.4377267822337301635744652424953398116882149670716141232569729223172939742940416733395e-1|\
a[13,2]=0.|\
a[13,3]=0.|\
a[13,4]=0.|\
a[13,5]=0.|\
a[13,6]=0.|\
a[13,7]=0.|\
a[13,8]=0.|\
a[13,9]=.6243650275201952087943586285809336252816312169030959172012504609444028241438248581173e-2|\
a[13,10]=.2000430971095773149944351654696478568290662322182649696087680691197048872391143823078|\
a[13,11]=-.8053283678049830368238571620489029119233928873370293148442058928084075077460302544840e-2|\
a[13,12]=.2115175280673965219157119035233996013168778251575505730512208770404786743066139905871e-1|\
a[14,1]=.2834992503635145630950235919207173122471376548964770977684956012393009143065795513785e-1|\
a[14,2]=0.|\
a[14,3]=0.|\
a[14,4]=0.|\
a[14,5]=0.|\
a[14,6]=0.|\
a[14,7]=0.|\
a[14,8]=0.|\
a[14,9]=.2491632048558174075389491488059951494598846535854176800982219995075912885766744587193e-2|\
a[14,10]=.2301387878545931496383998463737427687720871226381422342236583655735620108657836993957e-1|\
a[14,11]=-.3221559566929770987244760924671208781894636047606204610433085107190031098987004938258e-2|\
a[14,12]=.9884425494476646689463354144878852560408199827860146481292993078049373245839618405001e-2|\
a[14,13]=-.2130107713288873513843076428759273848866345654295724666320922464722154754985568313136e-1|\
a[15,1]=.3435118942902430010494322347351479430833531749807014262686507474123120416010457867571|\
a[15,2]=0.|\
a[15,3]=0.|\
a[15,4]=0.|\
a[15,5]=0.|\
a[15,6]=0.|\
a[15,7]=0.|\
a[15,8]=0.|\
a[15,9]=.2104519120236273856090970119990106557888074052256267000419050051487632641518018732685|\
a[15,10]=1.034274520572304119364829268288257099386679996983247401666929134177931632176349026735|\
a[15,11]=.6003036458644224870512404482066405749390780924061569454673075686417142117164254262878e-2|\
a[15,12]=.8559381250996195375780121060024077289150626526164160058172684354881277648341960563008|\
a[15,13]=-.9772350050367668108722648523725256330131076568928396776974412446349105799705851506077|\
a[15,14]=-.6600269804792946946162250138563276937205739812199748747775581736879654453322759683463|\
a[16,1]=-.1435740016721680695382063999350763666577559543783998809757153672896315044426183882232e-1|\
a[16,2]=0.|\
a[16,3]=0.|\
a[16,4]=0.|\
a[16,5]=0.|\
a[16,6]=0.|\
a[16,7]=0.|\
a[16,8]=0.|\
a[16,9]=-.3662532700490399702936857968489747917331190817335522078657913621382824038988807796287e-1|\
a[16,10]=.3502549756362136819768494069798465243467890824711035742020654749717518291597210559354e-1|\
a[16,11]=.3609460163621135089317866587583352398236899298642376718895880083960486970547825683491e-1|\
a[16,12]=-.2652199675536811063515959468346019236496270124574642848667252606942739787160130682669e-1|\
a[16,13]=.4456990113056981196389115375088399081043363230822267716707629092315111479614958673268e-1|\
a[16,14]=.1243430933313582432862255957417864480389734088951067419167759990001419776217292554191|\
a[16,15]=.4138296932394806944035124962043359604261929086744760344472227418812310333088685698274e-2|\
a[17,1]=.3560324044251202909756091163980891762641062223797488026536968101501275805051289760823|\
a[17,2]=0.|\
a[17,3]=0.|\
a[17,4]=0.|\
a[17,5]=0.|\
a[17,6]=0.|\
a[17,7]=0.|\
a[17,8]=0.|\
a[17,9]=-.450192758947562595966821779075956175110645100214763601190349|\
a[17,10]=.430527907083710898626656292808782917793030154094709462877146|\
a[17,11]=.5119730290110222376685569603940716920771257870306513863906805244405042813755411512463|\
a[17,12]=.9083036388864042603901591246381102139974962148199046305445452541870528153933236088278|\
a[17,13]=-1.239210933719339317573724691515340288544138892486057261860887966510000755220957594942|\
a[17,14]=-.6490486616717614651416723488790625539054028319671910976544025456235491510559878435372|\
a[17,15]=.2517089045868192922104805299489705414048878529314474912189256354259853776829630937658|\
a[17,16]=.7799064703455863988107567952823344760235405934115501870206452879298798513199886085571|\
a[18,1]=.1309356874065130664068812064188349801274704382131924878449566575565302965696195341197e-1|\
a[18,2]=0.|\
a[18,3]=0.|\
a[18,4]=0.|\
a[18,5]=0.|\
a[18,6]=0.|\
a[18,7]=0.|\
a[18,8]=0.|\
a[18,9]=0.|\
a[18,10]=0.|\
a[18,11]=0.|\
a[18,12]=0.|\
a[18,13]=-.9320530679851139459084619627671082378586315096846671421247697017556505173897578610165e-4|\
a[18,14]=.5053743342622993596400904431385907267709423447161223817027456630856526555478831396014e-1|\
a[18,15]=.8044703419444879791095791096101977976413118689308653610493721999399129417586629251430e-6|\
a[18,16]=.5917260294941711905287557427777172598443409719243215281782302034071342229921661278343e-3|\
a[18,17]=-.4016147221545573370646916849063755877322642479500938046774565993013424294867398455789e-6|\
a[19,1]=.2079264844660530125419445440007656521672552061443734079797586969853055549175505457737e-1|\
a[19,2]=0.|\
a[19,3]=0.|\
a[19,4]=0.|\
a[19,5]=0.|\
a[19,6]=0.|\
a[19,7]=0.|\
a[19,8]=0.|\
a[19,9]=0.|\
a[19,10]=0.|\
a[19,11]=0.|\
a[19,12]=0.|\
a[19,13]=.5826959188000859151019026978372841089514061030298715701031065480360641416298102920851e-3|\
a[19,14]=-.8017007323588159390833421865258527466405584659196335246554992680506588169863285718822e-2|\
a[19,15]=.4038476438471369403751708217435605704841172903308955066191655368223862388605213690921e-5|\
a[19,16]=.8546099980555061442250561145675356025101146220336224918025961310211940592009621595606e-1|\
a[19,17]=-.2044864809358042427067075696910043079044428375526774562331430989116458814609927891477e-5|\
a[19,18]=.1053285788244318933997994029790939973542409042351728431465827473723673651882417656762|\
a[20,1]=1.401534497957360214154462473557713067184864529175977331289881318884096354294079099114|\
a[20,2]=0.|\
a[20,3]=0.|\
a[20,4]=0.|\
a[20,5]=0.|\
a[20,6]=0.|\
a[20,7]=0.|\
a[20,8]=0.|\
a[20,9]=0.|\
a[20,10]=0.|\
a[20,11]=0.|\
a[20,12]=0.|\
a[20,13]=-.2302520009842212616162724103674156212611302982744556219175010157057031125814669239016|\
a[20,14]=-7.211068404669129056595822371068742471658564935099615697324849532576890894506619405031|\
a[20,15]=.3729015606948363352369953278521323402177595666786623882373057096229137360164435411243e-2|\
a[20,16]=-4.714154957271250206787781793922247570113233732218200980194845522013711035054762664884|\
a[20,17]=-.1763676575453492420538419950327976735749038866956001340593194717236122233799126229446e-2|\
a[20,18]=7.641305480386987655630293108802376511851733678139370059818519661401442202665741111270|\
a[20,19]=3.506020436597518349898960829497447109682129498933753736341591881470708008233521976557|\
a[21,1]=11.95146506941206867993723858307164016744736108265535168242754934626543968357331742096|\
a[21,2]=0.|\
a[21,3]=0.|\
a[21,4]=0.|\
a[21,5]=0.|\
a[21,6]=0.|\
a[21,7]=0.|\
a[21,8]=0.|\
a[21,9]=0.|\
a[21,10]=0.|\
a[21,11]=0.|\
a[21,12]=0.|\
a[21,13]=7.794809321081759687835167002317643882202842795989809549197917776161588225206322580459|\
a[21,14]=-56.45013938673257925235609911209042814404681000613405538635967763011214022629172907669|\
a[21,15]=.9123763069306449013445304492902766457096074504036737047499704936582270274950128398912e-1|\
a[21,16]=-12.73362799254348862019455243091992750381627175299189605168457824373779389828110581300|\
a[21,17]=-.3968959219047197123135428109397366747123830704331478729319411886202118671113516172493e-1|\
a[21,18]=54.43921418835708869962257651553077918614383784233053341001985423053366890118247056463|\
a[21,19]=-3.644116379215692368464069903613506458067214784092667356589342345057374050114156075061|\
a[21,20]=-.8045032499105099108990307879585794993156949132107878807481027183961246894903442258757|\
a[22,1]=-148.8094265071004884278388682686476255619306120821485965777899951377767737092911763254|\
a[22,2]=0.|\
a[22,3]=0.|\
a[22,4]=0.|\
a[22,5]=0.|\
a[22,6]=0.|\
a[22,7]=0.|\
a[22,8]=0.|\
a[22,9]=0.|\
a[22,10]=0.|\
a[22,11]=0.|\
a[22,12]=0.|\
a[22,13]=-91.72952782912564843579356624023216234952287290363542836291360346578688265538801398361|\
a[22,14]=707.6561449715983598345757192863357161548211289666495623584804744987957677893379157809|\
a[22,15]=-1.10563611857482440905296961311590930801338308942637769555540|\
a[22,16]=176.1345918838113725878598980760556604069995167623016865882869129962911416096097878945|\
a[22,17]=.4913848242148806622688983451644545574168846314027647925019604519368994965045299923826|\
a[22,18]=-684.2780004498149443582375356108950819560771678936002751371799726829821841834791232605|\
a[22,19]=27.99106049983982589842243321243804074460025184006686868209688958109916979926727384229|\
a[22,20]=13.19397100302823334436709643711532384350641596237449753683872220663989495376087330358|\
a[22,21]=1.251287812839804454501149741480560063172688300773964063605141347518040989702499199856|\
a[23,1]=-9.673079469481967636441261184332193958399514085718772596349277868068021458303626779169|\
a[23,2]=0.|\
a[23,3]=0.|\
a[23,4]=0.|\
a[23,5]=0.|\
a[23,6]=0.|\
a[23,7]=0.|\
a[23,8]=0.|\
a[23,9]=0.|\
a[23,10]=0.|\
a[23,11]=0.|\
a[23,12]=0.|\
a[23,13]=-4.469901508585055314438462277019603604978306814087514357488023393670679083633020106516|\
a[23,14]=45.51271286909526819682419504000527511789059078173984816890412459840121969200961260987|\
a[23,15]=-.713085086183826912791492024438246129930559805352394367050813e-1|\
a[23,16]=11.22736140684127415825906244799393842078268007767944830815221105133516977144595052189|\
a[23,17]=.1262443767176227245162379129091388093617868898191054263714925416869147773104813482457|\
a[23,18]=-43.54393395494833136058106249072421076238143044676214056937881652359375369765457150165|\
a[23,19]=.7871743075430589783987929949965509020645460914432340378113766124779028133099797867162|\
a[23,20]=.5322646967446842156693007086038866907853957768215038536520118921656033723449302296244|\
a[23,21]=.4224227339963253260102251274713887725750865388096033468497941673910509540050957057177|\
a[23,22]=.8591312495030671073084380314998594434411150562941549563989586466154235621165245563192e-1|\
a[24,1]=-10.06640324470547024033966069004268914722028247579687652710623604380152449409080444899|\
a[24,2]=0.|\
a[24,3]=0.|\
a[24,4]=0.|\
a[24,5]=0.|\
a[24,6]=0.|\
a[24,7]=0.|\
a[24,8]=0.|\
a[24,9]=-.3662532700490399702936857968489747917331190817335522078657913621382824038988807796287e-1|\
a[24,10]=.3502549756362136819768494069798465243467890824711035742020654749717518291597210559354e-1|\
a[24,11]=.3609460163621135089317866587583352398236899298642376718895880083960486970547825683491e-1|\
a[24,12]=-.2652199675536811063515959468346019236496270124574642848667252606942739787160130682669e-1|\
a[24,13]=-6.270889721814641435905531494788716038393561229573960230194057818533161624674313994502|\
a[24,14]=48.20792374425629890907021030081950639234925931416361161278899187780407980462426656808|\
a[24,15]=-.694471689136165640882395180583732834557754169149088630301342e-1|\
a[24,16]=12.68106902048502956983413709136098070661084838114121251454273060707937017246509534894|\
a[24,17]=.119671168968323754838161435501011294100927813964199613229864e-1|\
a[24,18]=-46.72497649924824080033582682426626955932013216597956070401309263301039263373634230581|\
a[24,19]=1.330296133266267113147100392982165913990335111912271192356479099067512051132965697343|\
a[24,20]=1.007667875033982983534389036199266577711627177936617199056121787956529680139072027935|\
a[24,21]=.2095120519336650916641223884754807028927707538644872411247284065032940106679251005781e-1|\
a[24,22]=.2101347063312641773177354243313964074244121884437574908902263894855162847478911411134e-1|\
a[24,23]=.9521960144171217941751015424545759073763602336583562405468424451848266905185171865534e-2|\
a[25,1]=-409.4780816777437087725890974093703576244243416067520683455326035855162023776088699896|\
a[25,2]=0.|\
a[25,3]=0.|\
a[25,4]=0.|\
a[25,5]=0.|\
a[25,6]=0.|\
a[25,7]=0.|\
a[25,8]=0.|\
a[25,9]=.2104519120236273856090970119990106557888074052256267000419050051487632641518018732685|\
a[25,10]=1.034274520572304119364829268288257099386679996983247401666929134177931632176349026735|\
a[25,11]=.6003036458644224870512404482066405749390780924061569454673075686417142117164254262878e-2|\
a[25,12]=.8559381250996195375780121060024077289150626526164160058172684354881277648341960563008|\
a[25,13]=-250.5169985474478604927776577293161303865840504207820779326393997812026874735614210230|\
a[25,14]=1946.424666523884277660537503282647585958298508957614274560610260899186136259514015246|\
a[25,15]=-3.045038821023103655061058090868608827869505440976021016842196622317831446605499698935|\
a[25,16]=490.6263795282817135212082652991680838415985422740616633051003594128766152337185220086|\
a[25,17]=1.566475895312709071154840670135974457395956152459667753199388690841173424714434871921|\
a[25,18]=-1881.974289940111733622172673770358706192159066384530557689275696031792911993357071098|\
a[25,19]=75.25922247248471752788377136433031498216206189142459440229301807516615379972994062700|\
a[25,20]=34.57343569803310676224343447365546896967286447935510158001529990937243976348724448442|\
a[25,21]=3.211476794409689614354173618470737551690229667488916278855754113243135684398993410117|\
a[25,22]=-.4604080417384143913072014042370588488672450952653828208427296561415079214017074427602|\
a[25,23]=-.8707183398418105224318841379579862457242520473889365722145748143125162133630944128398e-1|\
a[25,24]=-7.393518141583030675670169521955210639991857732491329543926346613193825315394087286297|\
a[26,1]=3.433474758535508789210934962575967811206238910720084588712755786644583035514752699598|\
a[26,2]=0.|\
a[26,3]=0.|\
a[26,4]=0.|\
a[26,5]=0.|\
a[26,6]=0.|\
a[26,7]=0.|\
a[26,8]=0.|\
a[26,9]=.2491632048558174075389491488059951494598846535854176800982219995075912885766744587193e-2|\
a[26,10]=.2301387878545931496383998463737427687720871226381422342236583655735620108657836993957e-1|\
a[26,11]=-.3221559566929770987244760924671208781894636047606204610433085107190031098987004938258e-2|\
a[26,12]=.9884425494476646689463354144878852560408199827860146481292993078049373245839618405001e-2|\
a[26,13]=2.162527993779225077883078419047573540457592253357327094851479956564246957314476133478|\
a[26,14]=-16.26998645464574213280656406601394890069875520402288517985775075363232756881970486667|\
a[26,15]=-.1285345021205245528435834174709350105380290375426545062302651848844352856037884822181|\
a[26,16]=-8.98915042666504253089307820833379330486511746063552853023189|\
a[26,17]=-.3485953632320253333870802018510136501924017672505137649688730136175086767654181319387e-2|\
a[26,18]=15.79361941133398075362351873886955741358533870251397376656158275266140525531011608606|\
a[26,19]=-.5744033309140950656281654820173358201483836631956754708231458398423255984252281047127|\
a[26,20]=-.3456020390213932966927224966081249825352372288276553067081833889419898565070467534157|\
a[26,21]=-.6622414902065850917316199913837577811330679927074186873906450413385445874036001388495e-2|\
a[26,22]=-.7777881292422041640325464586073643097593472096267591120155367761150273183248441708392e-2|\
a[26,23]=-.3560841924022749133388272326974373646752408187917065879526063406092336300493607300593e-2|\
a[26,24]=4.792825064499307996497977496298401894572969341393590555417712618624354747222657791607|\
a[26,25]=.153725464873068577844576387402512082757034273069877432944621|\
a[27,1]=32.30385208719854423269947344400315350913649750477846297617061421719281146058139852238|\
a[27,2]=0.|\
a[27,3]=0.|\
a[27,4]=0.|\
a[27,5]=0.|\
a[27,6]=-.317987696266205093901912847692712407988609169703103952205634e-2|\
a[27,7]=.8063977149061920772608217115203795063935431115674197501197468839656405367779525213500|\
a[27,8]=.9759831264123889790935228506842888513146720480030545503571875185550549213299958241991e-1|\
a[27,9]=.7785755781583989090275124464529272389997634605941819649588520345133050850477185489203|\
a[27,10]=.2048904238315994281894992020981056033120292350814206535748293420400885242747823516625|\
a[27,11]=-1.562615796274681883070709439505278252114628922364243608928053762634922556160297217820|\
a[27,12]=0.|\
a[27,13]=16.34298918823105706485042439739271747087533535041545512917666902744198799725970841669|\
a[27,14]=-154.5445552935436212307301896314710363993166836696091165017078152549564923882084122674|\
a[27,15]=1.569710887033348726920342834176217614662635935824970859658624964687079589089479471888|\
a[27,16]=3.276855450872481313214298172699007311655224049747336000450385269517693130775985884604|\
a[27,17]=-.5034892451936531763480407271997836265340810956916323972462042700071863164675818955838e-1|\
a[27,18]=153.3211518580416650705937678859146940112243631025945564907021486707139114294996134941|\
a[27,19]=7.175681863277204958467664848147841435678263080348653386540185145833155908488128910568|\
a[27,20]=-2.940367486753004819459176598969309892153205943807775979427615740476908865098135595635|\
a[27,21]=-.6658459460768031444707496760226288702819204931972568878708744783028558369468497032253e-1|\
a[27,22]=-.4623460549908436612292486685622172611769665140168592842374268449140643068786760618896e-1|\
a[27,23]=-.2041987335856794015393882286172697788485797748215817776751235910664984352284968100100e-1|\
a[27,24]=-53.35231064387358505159534411659981079740450904957915977996876390672711239156977103431|\
a[27,25]=-1.355487147150786549787321867059964040175545016141913251148206738329360142936656282958|\
a[27,26]=-1.571962758012327518829017351714592491776872191144425834618663282570958684038698495739|\
a[28,1]=-16.64514674863415128720312944039317587645603711308189782044257016154825923946758475845|\
a[28,2]=0.|\
a[28,3]=0.|\
a[28,4]=0.|\
a[28,5]=0.|\
a[28,6]=.5922327803245033080429900057980465247383895604442571368349896773084347972825775455007e-2|\
a[28,7]=.4703261599638411122172243032058941134553625307461088250108483236601604516650193568134|\
a[28,8]=.2996888638486790008539818370961923991368311216717812791841936858888827504094204242461|\
a[28,9]=-.2476568775939949146899922763298108258539580692639470955481886317480090967647905771626|\
a[28,10]=.1108950297714376828939998518390617145224451736006787182086245987785252503880550245038|\
a[28,11]=0.|\
a[28,12]=-.4917190438462291470706666287041940976780819072106730449888664749836403474888832394921|\
a[28,13]=-11.47431544272894969683894925643525363508424541308531757856483965863898534849416840511|\
a[28,14]=80.25931665762302725417024858864844001527933666235899875893849400507278534931158408231|\
a[28,15]=-.3841323039800428476253125267590291037469268413420882192068133107492120348263618466046|\
a[28,16]=7.281476674681075834713269509261361157676125818628777243483988994104498714011047355205|\
a[28,17]=-.1326993846122483795105717081760352748368273416167518843018178653526280269065470590467|\
a[28,18]=-81.07998325257307266746792897522552400060707166336329885641562357237166810196760593013|\
a[28,19]=-1.250374928356206395217681856561791199622537474924031863192434629401819729868852090550|\
a[28,20]=2.592635949695436810237763795043773249942264473592968880837586883560068434349818491911|\
a[28,21]=-.3014402983464045398301639972605268752644315372756414953420797074457552586137488110716|\
a[28,22]=.2213844607898323374517064515727737916952468390573184143179573617704323166985265217363|\
a[28,23]=.8275772747718929319559898709746931529962764354298098905497078729734353980896315305691e-1|\
a[28,24]=18.99606620406115204646724500372432639981751614122371589366718674999943569769696943522|\
a[28,25]=.2692319464096396856234680151283341674600519103489128451211866688910668614577677735665|\
a[28,26]=1.626748274470665374629893649296289339881250292841836802790201430504847697803528636395|\
a[28,27]=.4917190438462291470706666287041940976780819072106730449888664749836403474888832394921|\
a[29,1]=.8384798124090526646169687913728140859805331392249111310693346670107922625197375034871e-1|\
a[29,2]=0.|\
a[29,3]=0.|\
a[29,4]=0.|\
a[29,5]=0.|\
a[29,6]=-.1179493671009738143197550560312957753679619605907361507776128268875265788248790903515e-1|\
a[29,7]=-.2472990205688126523394738387431945983259928403533401326974984247503501083158412965835|\
a[29,8]=.9780808583677290122593130140812916655037406554767339407565991037499621093437371932341e-1|\
a[29,9]=.2175906892434206313600086517678603183441681200247821768799893467069296630467914197921|\
a[29,10]=0.|\
a[29,11]=.1375856067633252248656596321967877466474472229750848659754400903987833771639575727867|\
a[29,12]=.4398702297150466850587900923415450260461038902942613590425808839943205635447284745074e-1|\
a[29,13]=0.|\
a[29,14]=-.5137008137681933419570044566186303037387573636419640300869712169933398305905931343468|\
a[29,15]=.8263556911513155086442113083991534587014231586161685769224194977471882335420141183213|\
a[29,16]=25.70181397198118326258738829725199395111365563419600781824702737091645129169813134401|\
a[29,17]=0.|\
a[29,18]=0.|\
a[29,19]=0.|\
a[29,20]=0.|\
a[29,21]=0.|\
a[29,22]=0.|\
a[29,23]=0.|\
a[29,24]=-25.70181397198118326258738829725199395111365563419600781824702737091645129169813134401|\
a[29,25]=-.8263556911513155086442113083991534587014231586161685769224194977471882335420141183213|\
a[29,26]=.5137008137681933419570044566186303037387573636419640300869712169933398305905931343468|\
a[29,27]=-.4398702297150466850587900923415450260461038902942613590425808839943205635447284745074e-1|\
a[29,28]=-.1375856067633252248656596321967877466474472229750848659754400903987833771639575727867|\
a[30,1]=.1243805266540944128815164208687993162684914663596714231632892354628068537117612942798|\
a[30,2]=0.|\
a[30,3]=0.|\
a[30,4]=0.|\
a[30,5]=.2261202821975843014222386629792029011967523207426331439651447460281196206643404356021|\
a[30,6]=.1378858876180808806076958370164778145309694174914933853635428709475288586061552782365e-1|\
a[30,7]=-.6722101339966844497493995074143058569500863415253821828561997825320849038679063596730e-1|\
a[30,8]=0.|\
a[30,9]=0.|\
a[30,10]=-.8562389750854283547553497698795017721121215974115638028550665385850612741040225222977|\
a[30,11]=-1.963375228668589089282628500280938139881804405182674045535756631526916950083353845169|\
a[30,12]=-.2323328227241194012372462573089218472501081992304199949782180319905262045718872259601|\
a[30,13]=0.|\
a[30,14]=4.306607190864533494616689368765629477724325620534780926267640393608500758570100495873|\
a[30,15]=-2.927229632494654826597879112023904466876873949506336126307786635262992367484998786517|\
a[30,16]=-82.31316663978589444544923341054587077357619664281386893950601309356417181948645997040|\
a[30,17]=0.|\
a[30,18]=0.|\
a[30,19]=0.|\
a[30,20]=0.|\
a[30,21]=0.|\
a[30,22]=0.|\
a[30,23]=0.|\
a[30,24]=82.31316663978589444544923341054587077357619664281386893950601309356417181948645997040|\
a[30,25]=2.927229632494654826597879112023904466876873949506336126307786635262992367484998786517|\
a[30,26]=-4.306607190864533494616689368765629477724325620534780926267640393608500758570100495873|\
a[30,27]=.2323328227241194012372462573089218472501081992304199949782180319905262045718872259601|\
a[30,28]=1.963375228668589089282628500280938139881804405182674045535756631526916950083353845169|\
a[30,29]=.8562389750854283547553497698795017721121215974115638028550665385850612741040225222977|\
a[31,1]=.1034845616366797766729935465119103444997447982019713166066629728281981965079290745983|\
a[31,2]=0.|\
a[31,3]=0.|\
a[31,4]=.1220688873064072225896440828689620771395927148341621347412746563709055937325311521675|\
a[31,5]=.4825744903312466224751347801256881128659190238501680496794015023696413273862321544150|\
a[31,6]=-.3814096000156069997308862400056202056641130724784114774219699240039767479629669855696e-1|\
a[31,7]=0.|\
a[31,8]=-.5504995253108023241383885070205081774114143110000375617128363206424473498745141065969|\
a[31,9]=0.|\
a[31,10]=-.7119158115851892278876482620437943875782918824067455704957652139710574799878630163853|\
a[31,11]=-.5841296056715513404329887301584808720953353296452275957070524410065417676683463009109|\
a[31,12]=0.|\
a[31,13]=0.|\
a[31,14]=2.110463081258649321287173000466227503003750542789369878507182287710881470618943318741|\
a[31,15]=-.8374947367395721355257420230010379926952601753351235177405529298334532793741463162845e-1|\
a[31,16]=5.100214990723209140752959690433441131075450608628042491597346388445135412965217165555|\
a[31,17]=0.|\
a[31,18]=0.|\
a[31,19]=0.|\
a[31,20]=0.|\
a[31,21]=0.|\
a[31,22]=0.|\
a[31,23]=0.|\
a[31,24]=-5.100214990723209140752959690433441131075450608628042491597346388445135412965217165555|\
a[31,25]=.8374947367395721355257420230010379926952601753351235177405529298334532793741463162845e-1|\
a[31,26]=-2.110463081258649321287173000466227503003750542789369878507182287710881470618943318741|\
a[31,27]=0.|\
a[31,28]=.5841296056715513404329887301584808720953353296452275957070524410065417676683463009109|\
a[31,29]=.7119158115851892278876482620437943875782918824067455704957652139710574799878630163853|\
a[31,30]=.5504995253108023241383885070205081774114143110000375617128363206424473498745141065969|\
a[32,1]=.1933333333333333333333333333333333333333333333333333333333333333333333333333333333333|\
a[32,2]=0.|\
a[32,3]=.22|\
a[32,4]=-.8e-1|\
a[32,5]=0.|\
a[32,6]=0.|\
a[32,7]=.1099934255807247039194624048650683408451190582958464264636524271459687549994002654752|\
a[32,8]=-.2542970480762701613840685069971531221418356269767039208462421656164179875269042982442|\
a[32,9]=0.|\
a[32,10]=.8655707771166942543437703438210982818328474012330118593467368132762510892051242759318|\
a[32,11]=3.324164491140930831067995527865720183368600929369864071601998386039920635781409865040|\
a[32,12]=0.|\
a[32,13]=0.|\
a[32,14]=-12.01022233159779338823523851486618412603019426339968151272769528462035002110216728101|\
a[32,15]=.4766014662424932394304427768620618996029637820035802094825720242694315551196576125507|\
a[32,16]=-29.02430112210363905258026232136540995962512213324709106915239870601916450708546744075|\
a[32,17]=0.|\
a[32,18]=0.|\
a[32,19]=0.|\
a[32,20]=0.|\
a[32,21]=0.|\
a[32,22]=0.|\
a[32,23]=0.|\
a[32,24]=29.02430112210363905258026232136540995962512213324709106915239870601916450708546744075|\
a[32,25]=-.4766014662424932394304427768620618996029637820035802094825720242694315551196576125507|\
a[32,26]=12.01022233159779338823523851486618412603019426339968151272769528462035002110216728101|\
a[32,27]=0.|\
a[32,28]=-3.324164491140930831067995527865720183368600929369864071601998386039920635781409865040|\
a[32,29]=-.8655707771166942543437703438210982818328474012330118593467368132762510892051242759318|\
a[32,30]=.2542970480762701613840685069971531221418356269767039208462421656164179875269042982442|\
a[32,31]=-.1099934255807247039194624048650683408451190582958464264636524271459687549994002654752|\
a[33,1]=-.8333333333333333333333333333333333333333333333333333333333333333333333333333333333333|\
a[33,2]=1.388888888888888888888888888888888888888888888888888888888888888888888888888888888889|\
a[33,3]=0.|\
a[33,4]=0.|\
a[33,5]=-.75|\
a[33,6]=0.|\
a[33,7]=-.4925295437180263044226820491140213202002146815806577847190740839644346370048749342561|\
a[33,8]=0.|\
a[33,9]=0.|\
a[33,10]=0.|\
a[33,11]=0.|\
a[33,12]=0.|\
a[33,13]=0.|\
a[33,14]=0.|\
a[33,15]=0.|\
a[33,16]=0.|\
a[33,17]=0.|\
a[33,18]=0.|\
a[33,19]=0.|\
a[33,20]=0.|\
a[33,21]=0.|\
a[33,22]=0.|\
a[33,23]=0.|\
a[33,24]=0.|\
a[33,25]=0.|\
a[33,26]=0.|\
a[33,27]=0.|\
a[33,28]=0.|\
a[33,29]=0.|\
a[33,30]=0.|\
a[33,31]=.4925295437180263044226820491140213202002146815806577847190740839644346370048749342561|\
a[33,32]=.75|\
a[34,1]=.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111|\
a[34,2]=0.|\
a[34,3]=-.2222222222222222222222222222222222222222222222222222222222222222222222222222222222222|\
a[34,4]=0.|\
a[34,5]=0.|\
a[34,6]=0.|\
a[34,7]=0.|\
a[34,8]=0.|\
a[34,9]=0.|\
a[34,10]=0.|\
a[34,11]=0.|\
a[34,12]=0.|\
a[34,13]=0.|\
a[34,14]=0.|\
a[34,15]=0.|\
a[34,16]=0.|\
a[34,17]=0.|\
a[34,18]=0.|\
a[34,19]=0.|\
a[34,20]=0.|\
a[34,21]=0.|\
a[34,22]=0.|\
a[34,23]=0.|\
a[34,24]=0.|\
a[34,25]=0.|\
a[34,26]=0.|\
a[34,27]=0.|\
a[34,28]=0.|\
a[34,29]=0.|\
a[34,30]=0.|\
a[34,31]=0.|\
a[34,32]=0.|\
a[34,33]=.2222222222222222222222222222222222222222222222222222222222222222222222222222222222222|\
a[35,1]=.2858351403889715587960888421638364148529275378945964668924322897553490152559792262023|\
a[35,2]=.2916666666666666666666666666666666666666666666666666666666666666666666666666666666667|\
a[35,3]=.21875|\
a[35,4]=0.|\
a[35,5]=.1640625|\
a[35,6]=0.|\
a[35,7]=.2181943549455566583271882415813521070932888243221879411415164327116967439531911272777|\
a[35,8]=.1803928984786977668636352219467754377196200536418492285624347210514163759703679527180|\
a[35,9]=0.|\
a[35,10]=.2057138394048450188591207551229295422775700949828089053939914789386228504942804843989|\
a[35,11]=.2427157915817702399702829279594465157627459713866705419485763522859549196625913978401|\
a[35,12]=.2464657808136293058336092911818914077992281038693057051370210135284213379790417930740|\
a[35,13]=-3.449919407908908249798341546016226620603704606149316442883265523381128452524989278943|\
a[35,14]=.2288755621600360817607290607384585842942203725527402184592948392511281334278617959957|\
a[35,15]=.2832905997021514153215274190567333359784365954938557898314048426595070708424182066065|\
a[35,16]=3.210851258377666409601314905442367870055573203322387098512984999880577120008173123283|\
a[35,17]=-.2235387773648456999202337562141625079641252300836740320899016275445898395177373582441|\
a[35,18]=-.7071211572044190735187272862074872121300912319552061607910521928571247612111795934106|\
a[35,19]=3.211233451502870804081747292028565008932600344430223743249588034157195885590228893622|\
a[35,20]=1.409543483096697660304144743011231757690459455735489863573218752821178310978199657967|\
a[35,21]=-.1513620534437426131216022767425181110909630262036760559494590353712667648924754181285|\
a[35,22]=.3723505745270142764547240802146199843971210282021482987373568243836683323798121465643|\
a[35,23]=.2529787464063613367221999077621412859157757281294143192610824780367182739421617243696|\
a[35,24]=-3.210851258377666409601314905442367870055573203322387098512984999880577120008173123283|\
a[35,25]=-.2832905997021514153215274190567333359784365954938557898314048426595070708424182066065|\
a[35,26]=-.2288755621600360817607290607384585842942203725527402184592948392511281334278617959957|\
a[35,27]=-.2464657808136293058336092911818914077992281038693057051370210135284213379790417930740|\
a[35,28]=-.2427157915817702399702829279594465157627459713866705419485763522859549196625913978401|\
a[35,29]=-.2057138394048450188591207551229295422775700949828089053939914789386228504942804843989|\
a[35,30]=-.1803928984786977668636352219467754377196200536418492285624347210514163759703679527180|\
a[35,31]=-.2181943549455566583271882415813521070932888243221879411415164327116967439531911272777|\
a[35,32]=-.1640625|\
a[35,33]=-.21875|\
a[35,34]=-.2916666666666666666666666666666666666666666666666666666666666666666666666666666666667|\
b[1]=.1785714285714285714285714285714285714285714285714285714285714285714285714285714285714e-1|\
b[2]=.5859375e-2|\
b[3]=.1171875e-1|\
b[4]=0.|\
b[5]=.17578125e-1|\
b[6]=0.|\
b[7]=.234375e-1|\
b[8]=.29296875e-1|\
b[9]=0.|\
b[10]=.3515625e-1|\
b[11]=.41015625e-1|\
b[12]=.46875e-1|\
b[13]=0.|\
b[14]=.52734375e-1|\
b[15]=.5859375e-1|\
b[16]=.64453125e-1|\
b[17]=0.|\
b[18]=.1053521135717530196914960328878781622276730830805238840416702908213176249782427570033|\
b[19]=.1705613462417521823821203385538740858875554878027908047375010369442754416180982144816|\
b[20]=.2062293973293519407835264857011048947419142862595424540779715293772640762608018856579|\
b[21]=.2062293973293519407835264857011048947419142862595424540779715293772640762608018856579|\
b[22]=.1705613462417521823821203385538740858875554878027908047375010369442754416180982144816|\
b[23]=.1053521135717530196914960328878781622276730830805238840416702908213176249782427570033|\
b[24]=-.64453125e-1|\
b[25]=-.5859375e-1|\
b[26]=-.52734375e-1|\
b[27]=-.46875e-1|\
b[28]=-.41015625e-1|\
b[29]=-.3515625e-1|\
b[30]=-.29296875e-1|\
b[31]=-.234375e-1|\
b[32]=-.17578125e-1|\
b[33]=-.1171875e-1|\
b[34]=-.5859375e-2|\
b[35]=.1785714285714285714285714285714285714285714285714285714285714285714285714285714285714e-1|\
'

print(CopyRightInfo) 
print('\n'*2)
from PeterStone import * 
PrintFixed(a)
print('\n'*2)
print('# The estimate of the local truncation error is  ( 1/1000 )  h ( f(t1,x1)-f(t33,x33) )') 
print("print('k1=%e, k33=%e, dif=%e' % (k1, k33, abs(k1-k33)))")
print("\n")

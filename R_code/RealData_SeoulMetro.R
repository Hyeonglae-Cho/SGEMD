# Load data

station_info = read.csv('master_station_info.csv', header = TRUE)
station_dist = read.csv('master_station_distance.csv', header = TRUE)

# Obtain data during January 2021

station_info_01 = station_info[str_detect(station_info$date, '2021-01'),]
station_info_01_monthly.tmp = aggregate(total ~ name_kor, station_info_01, sum)

station_info_01_monthly = left_join(station_info_01_monthly.tmp, data.frame(name_kor = station_info_01$name_kor[1:243], lon = station_info_01$lon[1:243], lat = station_info_01$lat[1:243]), by = c('name_kor'))

# Obtain the adjacency matrix

station_adjacency_matrix.tmp = data.frame(matrix(rep(0, 243^2), nrow=243, ncol=243) ,row.names = station_info_01_monthly$name_kor)
colnames(station_adjacency_matrix.tmp) = station_info_01_monthly$name_kor

station_adjacency_matrix.tmp['가락시장', c('문정', '송파', '수서', '경찰병원')] = 1
station_adjacency_matrix.tmp['가산디지털단지', c('철산', '남구로')] = 1
station_adjacency_matrix.tmp['강남', c('교대', '역삼')] = 1
station_adjacency_matrix.tmp['강남구청', c('학동', '청담')] = 1
station_adjacency_matrix.tmp['강동', c('길동', '둔촌동', '천호')] = 1
station_adjacency_matrix.tmp['강동구청', c('천호', '몽촌토성')] = 1
station_adjacency_matrix.tmp['강변', c('구의', '잠실나루')] = 1
station_adjacency_matrix.tmp['개롱', c('오금', '거여')] = 1
station_adjacency_matrix.tmp['개화산', c('방화', '김포공항')] = 1
station_adjacency_matrix.tmp['거여', c('개롱', '마천')] = 1
station_adjacency_matrix.tmp['건대입구', c('구의', '성수', '뚝섬유원지', '어린이대공원')] = 1
station_adjacency_matrix.tmp['경복궁', c('안국', '독립문')] = 1
station_adjacency_matrix.tmp['경찰병원', c('가락시장', '오금')] = 1
station_adjacency_matrix.tmp['고덕', c('명일', '상일동')] = 1
station_adjacency_matrix.tmp['고려대', c('안암', '월곡')] = 1
station_adjacency_matrix.tmp['고속터미널', c('내방', '반포', '신반포', '사평', '잠원', '교대')] = 1
station_adjacency_matrix.tmp['공덕', c('대흥', '효창공원앞', '마포', '애오개')] = 1
station_adjacency_matrix.tmp['공릉', c('하계', '태릉입구')] = 1
station_adjacency_matrix.tmp['광나루', c('아차산', '천호')] = 1
station_adjacency_matrix.tmp['광명사거리', c('천왕', '철산')] = 1
station_adjacency_matrix.tmp['광화문', c('종로3가', '서대문')] = 1
station_adjacency_matrix.tmp['광흥창', c('상수', '대흥')] = 1
station_adjacency_matrix.tmp['교대', c('고속터미널', '남부터미널', '서초', '강남')] = 1
station_adjacency_matrix.tmp['구로디지털단지', c('대림', '신대방')] = 1
station_adjacency_matrix.tmp['구산', c('응암', '연신내')] = 1
station_adjacency_matrix.tmp['구의', c('건대입구', '강변')] = 1
station_adjacency_matrix.tmp['구파발', c('지축', '연신내')] = 1
station_adjacency_matrix.tmp['군자', c('장한평', '아차산', '어린이대공원', '중곡')] = 1

station_adjacency_matrix.tmp['굴포천', c('부평구청', '삼산체육관')] = 1 
station_adjacency_matrix.tmp['굽은다리', c('길동', '명일')] = 1
station_adjacency_matrix.tmp['금호', c('옥수', '약수')] = 1
station_adjacency_matrix.tmp['길동', c('강동', '굽은다리')] = 1
station_adjacency_matrix.tmp['길음', c('성신여대입구', '미아사거리')] = 1
station_adjacency_matrix.tmp['김포공항', c('개화산', '송정', '개화', '공항시장')] = 1
station_adjacency_matrix.tmp['까치산', c('화곡', '신정', '신정네거리')] = 1
station_adjacency_matrix.tmp['까치울', c('부천종합운동장', '온수')] = 1
station_adjacency_matrix.tmp['낙성대', c('서울대입구', '사당')] = 1
station_adjacency_matrix.tmp['남구로', c('가산디지털단지', '대림')] = 1
station_adjacency_matrix.tmp['남부터미널', c('교대', '양재')] = 1
station_adjacency_matrix.tmp['남성', c('숭실대입구', '총신대입구')] = 1
station_adjacency_matrix.tmp['남태령', c('사당', '선바위')] = 1
station_adjacency_matrix.tmp['남한산성입구', c('산성', '단대오거리')] = 1
station_adjacency_matrix.tmp['내방', c('총신대입구', '고속터미널')] = 1
station_adjacency_matrix.tmp['노원', c('중계', '마들', '창동', '상계')] = 1
station_adjacency_matrix.tmp['녹번', c('불광', '홍제')] = 1
station_adjacency_matrix.tmp['녹사평', c('삼각지', '이태원')] = 1
station_adjacency_matrix.tmp['논현', c('반포', '학동')] = 1
station_adjacency_matrix.tmp['단대오거리', c('남한산성입구', '신흥')] = 1
station_adjacency_matrix.tmp['답십리', c('마장', '장한평')] = 1
station_adjacency_matrix.tmp['당고개', c('상계', '별내별가람')] = 1
station_adjacency_matrix.tmp['당산', c('선유도', '국회의사당', '합정', '영등포구청')] = 1
station_adjacency_matrix.tmp['대림', c('신도림', '구로디지털단지', '남구로', '신풍')] = 1
station_adjacency_matrix.tmp['대청', c('일원', '학여울')] = 1
station_adjacency_matrix.tmp['대치', c('도곡', '학여울')] = 1
station_adjacency_matrix.tmp['대흥', c('공덕', '광흥창')] = 1
station_adjacency_matrix.tmp['도곡', c('매봉', '대치')] = 1
station_adjacency_matrix.tmp['도림천', c('신도림', '양천구청')] = 1
station_adjacency_matrix.tmp['도봉산', c('장암', '수락산')] = 1
station_adjacency_matrix.tmp['독립문', c('무악재', '경복궁')] = 1

station_adjacency_matrix.tmp['독바위', c('연신내', '불광')] = 1
station_adjacency_matrix.tmp['돌곶이', c('상월곡', '석계')] = 1
station_adjacency_matrix.tmp['동대문', c('혜화', '동대문역사문화공원', '종로5가', '동묘앞')] = 1
station_adjacency_matrix.tmp['동대문역사문화공원', c('동대문', '충무로', '을지로4가', '신당', '청구')] = 1
station_adjacency_matrix.tmp['동대입구', c('충무로', '약수')] = 1
station_adjacency_matrix.tmp['동묘앞', c('신설동', '동대문', '창신', '신당')] = 1
station_adjacency_matrix.tmp['동작', c('흑석', '신반포', '이촌', '총신대입구')] = 1
station_adjacency_matrix.tmp['둔촌동', c('강동', '올림픽공원')] = 1
station_adjacency_matrix.tmp['디지털미디어시티', c('증산', '월드컵경기장')] = 1
station_adjacency_matrix.tmp['뚝섬', c('한양대', '성수')] = 1
station_adjacency_matrix.tmp['뚝섬유원지', c('건대입구', '청담')] = 1
station_adjacency_matrix.tmp['마곡', c('송정', '발산')] = 1
station_adjacency_matrix.tmp['마들', c('수락산', '노원')] = 1
station_adjacency_matrix.tmp['마장', c('왕십리', '답십리')] = 1
station_adjacency_matrix.tmp['마천', c('거여')] = 1
station_adjacency_matrix.tmp['마포', c('여의나루', '공덕')] = 1
station_adjacency_matrix.tmp['마포구청', c('망원', '월드컵경기장')] = 1
station_adjacency_matrix.tmp['망원', c('마포구청', '합정')] = 1
station_adjacency_matrix.tmp['매봉', c('양재', '도곡')] = 1
station_adjacency_matrix.tmp['먹골', c('태릉입구', '중화')] = 1
station_adjacency_matrix.tmp['면목', c('상봉', '사가정')] = 1
station_adjacency_matrix.tmp['명동', c('회현', '충무로')] = 1
station_adjacency_matrix.tmp['명일', c('굽은다리', '고덕')] = 1
station_adjacency_matrix.tmp['모란', c('수진')] = 1
station_adjacency_matrix.tmp['목동', c('오목교', '신정')] = 1
station_adjacency_matrix.tmp['몽촌토성', c('강동구청', '잠실')] = 1
station_adjacency_matrix.tmp['무악재', c('홍제', '독립문')] = 1
station_adjacency_matrix.tmp['문래', c('영등포구청', '신도림')] = 1
station_adjacency_matrix.tmp['문정', c('가락시장', '장지')] = 1
station_adjacency_matrix.tmp['미사', c('하남풍산', '상일동')] = c(1, 0.5) #강일역이 없어서 바로 상일동역에 연결

station_adjacency_matrix.tmp['미아', c('미아사거리', '수유')] = 1
station_adjacency_matrix.tmp['미아사거리', c('미아', '길음')] = 1
station_adjacency_matrix.tmp['반포', c('고속터미널', '논현')] = 1
station_adjacency_matrix.tmp['발산', c('마곡', '우장산')] = 1
station_adjacency_matrix.tmp['방배', c('사당', '서초')] = 1
station_adjacency_matrix.tmp['방이', c('올림픽공원', '오금')] = 1
station_adjacency_matrix.tmp['방화', c('개화산')] = 1
station_adjacency_matrix.tmp['버티고개', c('한강진', '약수')] = 1
station_adjacency_matrix.tmp['보라매', c('신풍', '신대방삼거리')] = 1
station_adjacency_matrix.tmp['보문', c('안암', '창신')] = 1
station_adjacency_matrix.tmp['복정', c('산성', '장지')] = c(0.5, 1) #남위례역이 없어서 바로 산성 역으로 연결
station_adjacency_matrix.tmp['봉천', c('신림', '서울대입구')] = 1
station_adjacency_matrix.tmp['봉화산', c('신내', '화랑대')] = 1
station_adjacency_matrix.tmp['부천시청', c('상동', '신중동')] = 1
station_adjacency_matrix.tmp['부천종합운동장', c('까치울', '춘의')] = 1
station_adjacency_matrix.tmp['부평구청', c('산곡', '굴포천')] = 1
station_adjacency_matrix.tmp['불광', c('연신내', '녹번', '독바위', '역촌')] = 1
station_adjacency_matrix.tmp['사가정', c('용마산', '면목')] = 1
station_adjacency_matrix.tmp['사당', c('낙성대', '방배', '총신대입구', '남태령')] = 1
station_adjacency_matrix.tmp['산성', c('복정', '남한산성입구')] = c(0.5, 1) #남위례역이 없어서 바로 복정 역으로 연결
station_adjacency_matrix.tmp['삼각지', c('숙대입구', '녹사평', '신용산', '효창공원앞')] = 1
station_adjacency_matrix.tmp['삼산체육관', c('상동', '굴포천')] = 1
station_adjacency_matrix.tmp['삼성', c('선릉', '종합운동장')] = 1
station_adjacency_matrix.tmp['상계', c('노원', '당고개')] = 1
station_adjacency_matrix.tmp['상도', c('숭실대입구', '장승배기')] = 1
station_adjacency_matrix.tmp['상동', c('삼산체육관', '부천시청')] = 1
station_adjacency_matrix.tmp['상봉', c('면목', '중화')] = 1
station_adjacency_matrix.tmp['상수', c('합정', '광흥창')] = 1
station_adjacency_matrix.tmp['상왕십리', c('신당', '왕십리')] = 1
station_adjacency_matrix.tmp['상월곡', c('돌곶이', '월곡')] = 1
station_adjacency_matrix.tmp['상일동', c('고덕', '미사')] = c(1, 0.5) #강일역이 없어서 바로 미사역에 연결

station_adjacency_matrix.tmp['새절', c('응암', '증산')] = 1
station_adjacency_matrix.tmp['서대문', c('광화문', '충정로')] = 1
station_adjacency_matrix.tmp['서울대입구', c('봉천', '낙성대')] = 1
station_adjacency_matrix.tmp['서울역', c('회현', '숙대입구', '시청', '남영')] = 1
station_adjacency_matrix.tmp['서초', c('방배', '교대')] = 1
station_adjacency_matrix.tmp['석계', c('돌곶이', '태릉입구', '광운대', '신이문')] = 1
station_adjacency_matrix.tmp['석촌', c('잠실', '송파')] = 1
station_adjacency_matrix.tmp['선릉', c('삼성', '역삼')] = 1
station_adjacency_matrix.tmp['성수', c('뚝섬', '용답', '건대입구')] = 1
station_adjacency_matrix.tmp['성신여대입구', c('길음', '한성대입구')] = 1
station_adjacency_matrix.tmp['송정', c('김포공항', '마곡')] = 1
station_adjacency_matrix.tmp['송파', c('석촌', '가락시장')] = 1
station_adjacency_matrix.tmp['수락산', c('도봉산', '마들')] = 1
station_adjacency_matrix.tmp['수서', c('일원', '가락시장')] = 1
station_adjacency_matrix.tmp['수유', c('미아', '쌍문')] = 1
station_adjacency_matrix.tmp['수진', c('모란', '신흥')] = 1
station_adjacency_matrix.tmp['숙대입구', c('삼각지', '서울역')] = 1
station_adjacency_matrix.tmp['숭실대입구', c('상도', '남성')] = 1
station_adjacency_matrix.tmp['시청', c('충정로', '을지로입구', '종각', '서울역')] = 1
station_adjacency_matrix.tmp['신금호', c('청구', '행당')] = 1
station_adjacency_matrix.tmp['신길', c('영등포시장', '여의도', '영등포', '대방')] = 1
station_adjacency_matrix.tmp['신답', c('용두', '용답')] = 1
station_adjacency_matrix.tmp['신당', c('동대문역사문화공원', '상왕십리', '동묘앞', '청구')] = 1
station_adjacency_matrix.tmp['신대방', c('구로디지털단지', '신림')] = 1
station_adjacency_matrix.tmp['신대방삼거리', c('보라매', '장승배기')] = 1
station_adjacency_matrix.tmp['신도림', c('문래', '도림천', '대림')] = 1
station_adjacency_matrix.tmp['신림', c('신대방', '봉천')] = 1
station_adjacency_matrix.tmp['신사', c('잠원', '압구정')] = 1
station_adjacency_matrix.tmp['신설동', c('동묘앞', '제기동', '용두')] = 1
station_adjacency_matrix.tmp['신용산', c('삼각지', '이촌')] = 1
station_adjacency_matrix.tmp['신정', c('까치산', '목동')] = 1

station_adjacency_matrix.tmp['신정네거리', c('양천구청', '까치산')] = 1
station_adjacency_matrix.tmp['신중동', c('부천시청', '춘의')] = 1
station_adjacency_matrix.tmp['신촌', c('홍대입구', '이대')] = 1
station_adjacency_matrix.tmp['신풍', c('대림', '보라매')] = 1
station_adjacency_matrix.tmp['신흥', c('수진', '단대오거리')] = 1
station_adjacency_matrix.tmp['쌍문', c('수유', '창동')] = 1
station_adjacency_matrix.tmp['아차산', c('군자', '광나루')] = 1
station_adjacency_matrix.tmp['아현', c('이대', '충정로')] = 1
station_adjacency_matrix.tmp['안국', c('경복궁', '종로3가')] = 1
station_adjacency_matrix.tmp['안암', c('고려대', '보문')] = 1
station_adjacency_matrix.tmp['암사', c('천호')] = 1
station_adjacency_matrix.tmp['압구정', c('신사', '옥수')] = 1
station_adjacency_matrix.tmp['애오개', c('충정로', '공덕')] = 1
station_adjacency_matrix.tmp['약수', c('동대입구', '금호', '청구', '버티고개')] = 1
station_adjacency_matrix.tmp['양재', c('남부터미널', '매봉')] = 1
station_adjacency_matrix.tmp['양천구청', c('신정네거리', '도림천')] = 1
station_adjacency_matrix.tmp['양평', c('오목교', '영등포구청')] = 1
station_adjacency_matrix.tmp['어린이대공원', c('군자', '건대입구')] = 1
station_adjacency_matrix.tmp['여의나루', c('마포', '여의도')] = 1
station_adjacency_matrix.tmp['여의도', c('국회의사당', '샛강', '여의나루', '신길')] = 1
station_adjacency_matrix.tmp['역삼', c('강남', '선릉')] = 1
station_adjacency_matrix.tmp['역촌', c('응암', '불광')] = 1
station_adjacency_matrix.tmp['연신내', c('구파발', '불광', '구산', '독바위')] = 1
station_adjacency_matrix.tmp['영등포구청', c('양평', '영등포시장', '당산', '문래')] = 1
station_adjacency_matrix.tmp['영등포시장', c('영등포구청', '신길')] = 1
station_adjacency_matrix.tmp['오금', c('경찰병원', '방이', '개롱')] = 1
station_adjacency_matrix.tmp['오목교', c('목동', '양평')] = 1
station_adjacency_matrix.tmp['옥수', c('금호', '압구정')] = 1
station_adjacency_matrix.tmp['온수', c('까치울', '천왕', '역곡', '오류동')] = 1
station_adjacency_matrix.tmp['올림픽공원', c('둔촌동', '방이', '둔촌오륜', '한성백제')] = 1

station_adjacency_matrix.tmp['왕십리', c('상왕십리', '한양대', '마장', '행당')] = 1
station_adjacency_matrix.tmp['용답', c('신답', '성수')] = 1
station_adjacency_matrix.tmp['용두', c('신설동', '신답')] = 1
station_adjacency_matrix.tmp['용마산', c('사가정', '중곡')] = 1
station_adjacency_matrix.tmp['우장산', c('발산', '화곡')] = 1
station_adjacency_matrix.tmp['월곡', c('고려대', '상월곡')] = 1
station_adjacency_matrix.tmp['월드컵경기장', c('디지털미디어시티', '마포구청')] = 1
station_adjacency_matrix.tmp['을지로3가', c('을지로4가', '을지로입구', '종로3가', '충무로')] = 1
station_adjacency_matrix.tmp['을지로4가', c('을지로3가', '동대문역사문화공원', '종로3가')] = 1
station_adjacency_matrix.tmp['을지로입구', c('시청', '을지로3가')] = 1
station_adjacency_matrix.tmp['응암', c('역촌', '구산', '새절')] = 1
station_adjacency_matrix.tmp['이대', c('신촌', '아현')] = 1
station_adjacency_matrix.tmp['이촌', c('신용산', '동작')] = 1
station_adjacency_matrix.tmp['이태원', c('녹사평', '한강진')] = 1
station_adjacency_matrix.tmp['일원', c('대청', '수서')] = 1
station_adjacency_matrix.tmp['잠실', c('잠실새내', '잠실나루', '석촌', '몽촌토성')] = 1
station_adjacency_matrix.tmp['잠실나루', c('잠실', '강변')] = 1
station_adjacency_matrix.tmp['잠실새내', c('잠실', '종합운동장')] = 1
station_adjacency_matrix.tmp['잠원', c('신사', '고속터미널')] = 1
station_adjacency_matrix.tmp['장승배기', c('상도', '신대방삼거리')] = 1
station_adjacency_matrix.tmp['장암', c('도봉산')] = 1
station_adjacency_matrix.tmp['장지', c('문정', '복정')] = 1
station_adjacency_matrix.tmp['장한평', c('답십리', '군자')] = 1
station_adjacency_matrix.tmp['제기동', c('신설동', '청량리')] = 1
station_adjacency_matrix.tmp['종각', c('종로3가', '시청')] = 1
station_adjacency_matrix.tmp['종로3가', c('광화문', '을지로4가', '종각', '종로5가', '안국', '을지로3가')] = 1
station_adjacency_matrix.tmp['종로5가', c('종로3가', '동대문')] = 1
station_adjacency_matrix.tmp['종합운동장', c('잠실새내', '삼성', '봉은사', '삼전')] = 1
station_adjacency_matrix.tmp['중계', c('노원', '하계')] = 1

station_adjacency_matrix.tmp['중곡', c('용마산', '군자')] = 1
station_adjacency_matrix.tmp['중화', c('먹골', '상봉')] = 1
station_adjacency_matrix.tmp['증산', c('디지털미디어시티', '새절')] = 1
station_adjacency_matrix.tmp['지축', c('삼송', '구파발')] = 1
station_adjacency_matrix.tmp['창동', c('녹천', '방학', '쌍문', '노원')] = 1
station_adjacency_matrix.tmp['창신', c('보문', '동묘앞')] = 1
station_adjacency_matrix.tmp['천왕', c('온수', '광명사거리')] = 1
station_adjacency_matrix.tmp['천호', c('암사', '강동구청', '광나루', '강동')] = 1
station_adjacency_matrix.tmp['철산', c('광명사거리', '가산디지털단지')] = 1
station_adjacency_matrix.tmp['청구', c('신당', '약수', '신금호', '동대문역사문화공원')] = 1
station_adjacency_matrix.tmp['청담', c('뚝섬유원지', '강남구청')] = 1
station_adjacency_matrix.tmp['청량리', c('회기', '제기동')] = 1
station_adjacency_matrix.tmp['총신대입구', c('내방', '남성', '동작', '사당')] = 1
station_adjacency_matrix.tmp['춘의', c('신중동', '부천종합운동장')] = 1
station_adjacency_matrix.tmp['충무로', c('동대문역사문화공원', '명동', '을지로3가', '동대입구')] = 1
station_adjacency_matrix.tmp['충정로', c('서대문', '애오개', '아현', '시청')] = 1
station_adjacency_matrix.tmp['태릉입구', c('석계', '화랑대', '공릉', '먹골')] = 1
station_adjacency_matrix.tmp['하계', c('중계', '공릉')] = 1
station_adjacency_matrix.tmp['하남풍산', c('미사', '하남시청')] = 1
station_adjacency_matrix.tmp['학동', c('강남구청', '논현')] = 1
station_adjacency_matrix.tmp['학여울', c('대치', '대청')] = 1
station_adjacency_matrix.tmp['한강진', c('이태원', '버티고개')] = 1
station_adjacency_matrix.tmp['한성대입구', c('혜화', '성신여대입구')] = 1
station_adjacency_matrix.tmp['한양대', c('왕십리', '뚝섬')] = 1
station_adjacency_matrix.tmp['합정', c('홍대입구', '당산', '상수', '망원')] = 1
station_adjacency_matrix.tmp['행당', c('왕십리', '신금호')] = 1
station_adjacency_matrix.tmp['혜화', c('한성대입구', '동대문')] = 1
station_adjacency_matrix.tmp['홍대입구', c('신촌', '합정')] = 1
station_adjacency_matrix.tmp['홍제', c('녹번', '무악재')] = 1
station_adjacency_matrix.tmp['화곡', c('우장산', '까치산')] = 1
station_adjacency_matrix.tmp['화랑대', c('태릉입구', '봉화산')] = 1
station_adjacency_matrix.tmp['회현', c('명동', '서울역')] = 1
station_adjacency_matrix.tmp['효창공원앞', c('공덕', '삼각지')] = 1

station_adjacency_matrix.tmp = station_adjacency_matrix.tmp[,1:243]
rownames(station_adjacency_matrix.tmp) = NULL
colnames(station_adjacency_matrix.tmp) = NULL

station_adjacency_matrix = Matrix(as.matrix(station_adjacency_matrix.tmp), sparse=TRUE)

# Define graph

g_station = graph_from_adjacency_matrix(station_adjacency_matrix, mode=c('undirected'), weighted=TRUE)

graph_attr(g_station, 'xy') = station_info_01_monthly[,3:4] # Coordinates of 243 stations (vertices)
colnames(g_station$xy) = c('x', 'y')
graph_attr(g_station, 'sA') = station_adjacency_matrix

# Define graph signal

g_station_signal.standard = (station_info_01_monthly$total - mean(station_info_01_monthly$total)) / sd(station_info_01_monthly$total)
g_station_signal = g_station_signal.standard

# Graph EMD
g_station.GEMD = GEMD(g_station, g_station_signal, 4)

# Statistical graph EMD
g_station.SGEMD = GEMD.refl.GFT.ebayesthresh(g_station, g_station_signal, 4, 'KnnAvg', 0.02, 5, connection = 'neighbor')

# Define denoised part for each method

g_station.GEMD.denoised = g_station.GEMD[[2]] + g_station.GEMD[[3]] + g_station.GEMD[[4]] + g_station.GEMD[[5]]
g_station.SGEMD.denoised = g_station.SGEMD[[2]] + g_station.SGEMD[[3]] + g_station.SGEMD[[4]] + g_station.SGEMD[[5]]

# Define color range for visualization

color_range = c(g_station_signal, g_station.GEMD[[1]], g_station.GEMD.denoised, g_station.SGEMD[[1]], g_station.SGEMD.denoised)

# Plot graph signal

plot_signal1(g_station, color_range, g_station_signal, size = 3, custom_colours = custom_colours1)

# Visualization of the results from GEMD

plot_signal1(g_station, color_range, g_station.GEMD[[1]], size = 3, custom_colours = custom_colours1) # Extracted noise by GEMD
plot_signal1(g_station, c(-1, 1), g_station.GEMD[[2]], size = 3, custom_colours = custom_colours1) # 1st IMF decomposed by GEMD
plot_signal1(g_station, c(-1, 1), g_station.GEMD[[3]], size = 3, custom_colours = custom_colours1) # 2nd IMF decomposed by GEMD
plot_signal1(g_station, c(-1, 1), g_station.GEMD[[4]], size = 3, custom_colours = custom_colours1) # 3rd IMF decomposed by GEMD
plot_signal1(g_station, c(g_station.GEMD[[5]], g_station.SGEMD[[5]]), g_station.GEMD[[5]], size = 3, custom_colours = custom_colours1) # Residue decomposed by GEMD

plot_signal1(g_station, color_range, g_station.GEMD.denoised, size = 3, custom_colours = custom_colours1) # Denoised part decomposed by GEMD

# Visualization of the results from SGEMD
plot_signal1(g_station, color_range, g_station.SGEMD[[1]], size = 3, custom_colours = custom_colours1) # Extracted noise by SGEMD
plot_signal1(g_station, c(-1, 1), g_station.SGEMD[[2]], size = 3, custom_colours = custom_colours1) # 1st IMF decomposed by SGEMD
plot_signal1(g_station, c(-1, 1), g_station.SGEMD[[3]], size = 3, custom_colours = custom_colours1) # 2nd IMF decomposed by SGEMD
plot_signal1(g_station, c(-1, 1), g_station.SGEMD[[4]], size = 3, custom_colours = custom_colours1) # 3rd IMF decomposed by SGEMD
plot_signal1(g_station, c(g_station.GEMD[[5]], g_station.SGEMD[[5]]), g_station.SGEMD[[5]], size = 3, custom_colours = custom_colours1) # Residue decomposed by SGEMD

plot_signal1(g_station, color_range, g_station.SGEMD.denoised, size = 3, custom_colours = custom_colours1) # Denoised part decomposed by SGEMD

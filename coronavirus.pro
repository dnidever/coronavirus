pro coronavirus

;; Data from WHO coronavirus daily situation reports
;; https://www.who.int/emergencies/diseases/novel-coronavirus-2019/situation-reports
  
str = importascii('coronavirus.txt',/header)
usdaily = importascii('daily.csv',delim=',',/header)  ; from covidtracking.com
; add the covidtracking data to STR
add_tag,str,'us2',0L,str
; convert to MMDDYYYY
date2 = strmid(strtrim(usdaily.date,2),4,4)+strmid(strtrim(usdaily.date,2),0,4)
match,str.date,date2,ind1,ind2,/sort
str[ind1].us2 = long(usdaily[ind2].positiveincrease)

setdisp
!p.font = 0
file='coronavirus'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=10.0,ysize=8.5
xr = [1,max(str.num)+14]
yr = [1,max(str.all)*7]
plot,[1],/nodata,xr=xr,yr=yr,xs=1,ys=1,/ylog,xtit='Date',ytit='New Confirmed Cases',$
     charsize=1.5,title='New Confirmed Coronavirus Cases',$
     xminor=4,xticks=4,xtickv=[1,20,40,60,80],xtickn=['Jan 21','Feb 9','Feb 29','Mar 20','Apr 9']

;; ---- World minus China ----
gdw = where(str.all gt 0 and str.num ge 30,ngdw)
wcoef = robust_poly_fit(str[gdw].num,alog10(str[gdw].all),1)
g1 = where(str.all gt 0)
oplot,str[g1].num,str[g1].all,ps=8,co=250
x = findgen(100)+28
;oplot,x,10^poly(x,wcoef),co=250,thick=3 ;,linestyle=2
;; doubling time
wdouble = alog10(2)/wcoef[1]
;xyouts,5,6000,'Doubling Times:',align=0,charsize=1.7,co=0
;xyouts,5,3500,stringize(wdouble,ndec=1)+' days',align=0,charsize=1.7,co=250

; fit logistic curve to Italy data
initpar = [1e6,80.0,0.50]
g1b = where(str.all gt 0 and str.num gt 40)
fpar1 = mpfitfun('func_logisticderivlog',str[g1b].num,alog10(str[g1b].all),str[g1b].num*0+1,initpar)
x = findgen(100)+40
m = func_logisticderiv(x,fpar1)
;oplot,x,m,co=250,thick=5,linestyle=2


;; ---- Italy ----
green = fsc_color('forest green',1)
gdit = where(str.italy gt 0 and str.num gt 52,ngdit)
itcoef = robust_poly_fit(str[gdit].num,alog10(str[gdit].italy),1)
g2 = where(str.italy gt 0)
oplot,str[g2].num,str[g2].italy,ps=8,co=green,sym=1.0
x = findgen(100)+52
;oplot,x,10^poly(x,itcoef),co=green,thick=1 ;,linestyle=2
;; doubling time
itdouble = alog10(2)/itcoef[1]
;xyouts,5,1100,stringize(itdouble,ndec=1)+' days',align=0,charsize=1.7,co=green

; fit logistic curve to Italy data
;initpar = [150000.,68.0,0.15]
initpar = [200000.,69.0,0.20]
g2b = where(str.italy gt 0 and str.num gt 40)
;fpar2 = mpfitfun('func_logisticderiv',str[g2b].num,str[g2b].italy,str[g2b].num*0+1,initpar)
fpar2 = mpfitfun('func_logisticderivlog',str[g2b].num,alog10(str[g2b].italy),str[g2b].num*0+1,initpar)
x = findgen(100)+40
m = func_logisticderiv(x,fpar2)
;oplot,x,m,co=green,thick=5,linestyle=2


;; ---- US ----
gdus = where(str.us gt 0 and str.num ge 44,ngdus)
uscoef = robust_poly_fit(str[gdus].num,alog10(str[gdus].us),1)
g2 = where(str.us gt 0)
plotsym,0,1.2,thick=5
oplot,str[g2].num,str[g2].us,ps=8,co=70 ;,sym=1.2
x = findgen(100)+44
oplot,x,10^poly(x,uscoef),co=80,thick=3 ;,linestyle=2
;; covidtracking.com data
g3 = where(str.us2 gt 0)
oplot,str[g3].num,str[g3].us2,ps=1,co=50,sym=1.3,thick=5
;; doubling time
usdouble = alog10(2)/uscoef[1]
;xyouts,5,2000,stringize(usdouble,ndec=1)+' days (last 42 days)',align=0,charsize=1.7,co=50
oplot,[12],[4e4],ps=8,co=70
xyouts,13,3.6e4,'WHO/CDC',align=0,co=70,charsize=1.0
oplot,[26],[4e4],ps=1,co=70,sym=1.2
xyouts,27,3.6e4,'covidtracking.com',align=0,co=70,charsize=1.0

;; exponential, just the last week
g4 = where(str.num ge 65 and str.us2 gt 0,ng4)
uscoef_thisweek = robust_poly_fit(str[g4].num,alog10(str[g4].us2),1)
x = findgen(30)+65
;oplot,x,10^poly(x,uscoef_thisweek),co=90,thick=3
usdouble_thisweek = alog10(2)/uscoef_thisweek[1]
;xyouts,5,1100,stringize(usdouble_thisweek,ndec=1)+' days (last 22 days)',align=0,charsize=1.7,co=90

; logistic curve
initpar = [150000.,80.0,0.30]
;initpar = [4e4,80.0,fpar2[2]]
;fpar = mpfitfun('func_logisticderivlog',str[g3].num,alog10(str[g3].us2),str[g3].num*0+1,initpar)
g4 = where(str.num ge 58 and str.us2 gt 0,ng4)
;g4 = where(str.us2 gt 0,ng4)
parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},3)
parinfo[2].fixed = 1
fpar3 = mpfitfun('func_logisticderivlog',str[g4].num,alog10(str[g4].us2),str[g4].us2*0+1,initpar)
;fpar3 = mpfitfun('func_logisticderivlog',str[g4].num,alog10(str[g4].us2),sqrt(str[g4].us2)>1,initpar)
;fpar3 = mpfitfun('func_logisticderiv',str[g4].num,str[g4].us2,sqrt(str[g4].us2)>1,initpar,parinfo=parinfo)
;fpar3 = mpfitfun('func_logisticderivlog',str[g4].num,alog10(str[g4].us2),str[g4].num*0+1,initpar,parinfo=parinfo)
x = findgen(100)+58
m = func_logisticderiv(x,fpar3)
;oplot,x,m,co=80,thick=5,linestyle=2


;; lines for 1,000 cases
gray = fsc_color('gray',2)
oplot,[0,56.5],[1000,1000],linestyle=1
oplot,[56.5,56.5],[1,1000],linestyle=1
;oplot,[56.5,56.5],[5,1000],linestyle=1
;oplot,[56.5,56.5],[1,2.1],linestyle=1

;xyouts,56.5,4,'1 day to',align=0.5,co=0,charsize=1.1
;xyouts,56.5,3,'1000 new US',align=0.5,co=0,charsize=1.1
;xyouts,56.5,2.3,'cases a day',align=0.5,co=0,charsize=1.1

;; lines for 10,000 cases a new
g = first_el(where(10^poly(x,uscoef) ge 1e4,ng))
x1e4 = x[g[0]]
oplot,[0,x1e4],[10000,10000],linestyle=1
oplot,[0,0]+x1e4,[1,1e4],linestyle=1
;oplot,[0,0]+x1e4,[1,6],linestyle=1
;oplot,[0,0]+x1e4,[30,10000],linestyle=1

;ndays1e4 = x1e4-max(str.num)
;xyouts,63.8,20,strtrim(long(ndays1e4),2)+' days to',align=0.5,co=0,charsize=1.1
;xyouts,63.8,13,'10,000 new US',align=0.5,co=0,charsize=1.1
;xyouts,63.8,8.5,'cases a day',align=0.5,co=0,charsize=1.1

;; lines for 100,000 cases a new
g = first_el(where(10^poly(x,uscoef) ge 1e5,ng))
x1e5 = x[g[0]]
oplot,[0,100],[1e5,1e5],linestyle=1
;oplot,[0,x1e5],[1e5,1e5],linestyle=1
;oplot,[0,0]+x1e5,[1,60],linestyle=1
;oplot,[0,0]+x1e5,[300,1e5],linestyle=1

;ndays1e5 = x1e5-max(str.num)
;xyouts,72.5,200,strtrim(long(ndays1e5),2)+' days to',align=0.5,co=0,charsize=1.1
;xyouts,72.5,130,'100,000 new US',align=0.5,co=0,charsize=1.1
;xyouts,72.5,85,'cases a day',align=0.5,co=0,charsize=1.1

;; Today's date
jd = systime(/julian)
caldat,jd,month,day,year,hour,minute,second
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
datestr = months[month-1]+' '+strtrim(long(day),2)+', '+strtrim(long(year),2)
xyouts,5,yr[1]*0.45,datestr,align=0,charsize=1.8,co=0

;;xyouts,8.5,420,'Logistic',align=0,charsize=1.7
;;oplot,[5,8],[500,500],linestyle=2,co=0,thick=6
;xyouts,8.5,4200,'Logistic',align=0,charsize=1.7
;oplot,[5,8],[5000,5000],linestyle=2,co=0,thick=6

legend_old,['World minus China','US','Italy'],textcolor=[250,70,green],charsize=1.7,box=0,pos=[2,90000L]
;legend_old,[datestr,'','World minus China','US'],textcolor=[0,255,250,70],charsize=1.7,box=0,/top,/left
;pos=[1,8000]

ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell

stop

end

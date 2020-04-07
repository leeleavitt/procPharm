# Mac OS procPharm Setup

1.  Install [XQuartz-2.7.11.dmg](https://dl.bintray.com/xquartz/downloads/XQuartz-2.7.11.dmg/)
2.  Install [GTK_2.24.17-X11.pkg](http://r.research.att.com/libs/GTK_2.24.17-X11.pkg)
3.  Now Open `R` and type
````
install.packages('cairoDevice')
````
  * Make sure to click *Install Dependencies*

  * R should prompt you to select a CRAN location. I used USA (OR) just because it is close to Utah


### A couple of notes
Please do not go fullscreen (ie the green button) it causes the RDViewer and other functions to break and as of right now YOU’LL HAVE TO START FROM YOUR LAST SAVE
 \
It runs a little slower than on the computers at the lab, but be patient (especially with `ROIReview`)


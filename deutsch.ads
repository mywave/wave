!
!	@(#)deutsch.ads 04/01/06 14:21:20 1.4
! $Log: deutsch.ads,v $
! Revision 1.9  2005/03/04 09:57:43  frkochw
! Einbinden von neuronalen Netzen zur Windbestimmung.
!
! Revision 1.8  2005/02/02 15:12:46  frkochw
! Ocean Weather Windverzeichnis.
!
! Revision 1.7  2004/10/26 16:03:30  frkochw
! PolRatio eingefügt.
!
! Revision 1.6  2004/10/25 12:52:57  frkochw
! Voreinstellung für die Landseemaske ist jetzt lokal.
!
! Revision 1.5  2004/09/27 12:53:32  frkochw
! Copyright Notiz eingefügt.
!
! Revision 1.4  2004/02/26 08:39:39  frkochw
! Neue Items im Druckmenue von modell.pro: Scalloping, Farbkodierung und eindeutige Richtungen ncah Wackermann.
!
! Revision 1.3  2004/02/18 16:34:14  frkochw
! Zusätzliche Items für automatische Richtungsfindung.
!
! Revision 1.2  2004/01/06 14:03:48  frkochw
! Übergang auf CVS
!
!
Copyright:		GKSS SAR-Wind Program WiSAR, COPYRIGHT (C) 2004 BY GKSS\n\nTHIS PROGRAM IS EXPERIMENTAL AND IS PROVIDED "AS IS" WITHOUT\nREPRESENTATION OF WARRANTY OF ANY KIND, EITHER EXPRESS OR\nIMPLIED. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF\nTHE PROGRAM IS WITH THE USER.
Ausgabeformat:		Ausgabeformat: 
Hintergrund:		Hintergrund: 
Skalierung:		Skalierung: 
geglaettet:		geglättet
Balken:			Balken
Windgeschwindigkeit:	Windgeschwindigkeit
Keiner:			Keiner
KeineRichtungen:	Keine Richtungen gefunden
Mischen:		Mischen
Grau:			Grau
Spektrale:		Spektrale Richtungen
Scalloping:		Scalloping
Farbkodierung:		Farbkodierung
LaengenBreiten:		Längen und Breiten
FFTRichtungen:		FFT Richtungen
EindeutigeRichtungen:	Eindeutige Richtungen
Drucken:		Drucken
Analysenauswahl:	 Analysenauswahl
Clear:			Lösche 
Modell:			Modell
Heading:		Heading: 
ERSDaten:		/h/frkochw/ori/
OKohneLG:		OK ohne LG
BandBild:		 Bandnummer-Bildnummer <ENTER> eingeben 
KeinStreifenName:	Kein Name für den SAR Streifen gewählt
Ebene:			 Ebene
Palette:		/h/frkochw/paletten/rainbow.hdf
noSUNraster:		Dies ist kein SUN Rasterbild !
oirms:			/h/frkochw/02/vs/
Inhalt:			/h/frkochw/ERS/Band
Wind:			 Windhistogramm
WindVerteilung:		Windgeschwindigkeit über Wasser
Gegenwind:		Gegenwind 
Querwind:		Querwind 
anpat1:			/home/kochw/routine/openimop/anpat1.dat
anpat2:			/home/kochw/routine/openimop/anpat2.dat
anpat2b:		/home/kochw/routine/openimop/anpat2b.dat
anpat2c:		/home/kochw/routine/openimop/anpat2c.dat
anpat3:			/home/kochw/routine/openimop/anpat3.dat
anpat3b:		/home/kochw/routine/openimop/anpat3b.dat
anpat3c:		/home/kochw/routine/openimop/anpat3c.dat
anpat4:			/home/kochw/routine/openimop/anpat4.dat
anpat5:			/home/kochw/routine/openimop/anpat5.dat
powerloss1:		/home/kochw/routine/openimop/poloco1.dat
powerloss2:		/home/kochw/routine/openimop/poloco21.dat
RADARSATDaten:		/net/gmssun3/export/home/RADARSAT/
RADARSATOutput:		/net/gmssun3/export/home2/ScanSAR/ori
ASAR_XCA:		/home/kochw/xca/ASA_XCA_AX*
ENVISATOutput:		/net/gmssun16/export/home/frkochw/envisat/
Messungen:		/h/gmssun3-h2/frkochw/remo/messungen/messungen 
ReMoWind:		/h/gmssun3-h2/frkochw/remo/remo000xe
HiRLAMWind:		/h/gmssun3-h2/ScanSAR/HIRLAM_DAT/
HiRLAMWindEs:		/h/gmssun3-h2/frkochw/scansar/v
NOGABWind:		/h/gmssun3-h3/APL/
SCATWind:		/h/gmssun3-h3/APL/
DWDWind:		/net/gmssun16/export/home1/ENVISAT/DWD_DAT/
OWWind:			/h/frkochw/tmp/jh/
MELWind:		/home/kochw/nogab/
DWDLMWind:		/home/kochw/dwd/
Orbit:			Orbit: 
Richtungeindeutigmachen:Richtungen eindeutig machen (WK)
Wacker:			Richtungen eindeutig machen (Chris)
Druckausrechnen:	Bodendruck ausrechnen
Modellanzeigen:		Modell anzeigen
Tabelleausgeben:	Tabelle ausgeben
Datensichern:		Daten sichern
Hintergrundbild:	 Hintergrundbild
Rechnername:		echo local
MessungenOutput:	/net/gmssun21/export/home1/Wackerman/DAT/ori/
CBandModell:		smo5v
NN:			/home/kochw/openimop/NN_beta/16x14x12x10x8_51.0.net
NNprog:			/home/kochw/openimop/NN_beta/calc2v_NN
LandSeeDaten:		/home/kochw/gmt/tmp
LandSeeDatenErstellen:	rsh gmssun10 cd gmt;cdf.x
PolRatio:		model1

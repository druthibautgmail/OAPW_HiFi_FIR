# OAPW_HiFi_AL

## Open Ambiophonics Processing Wrapper – ALSA Loopback Edition

**OAPW_HiFi_AL** ist eine für den Raspberry Pi optimierte Echtzeit-Audiosignalverarbeitung in C++. Die Software implementiert den *Recursive Ambiophonic Crosstalk Elimination* (RACE)-Algorithmus, basierend auf den wegweisenden Arbeiten von Ralph Glasgal und dem Ambiophonics Institute.

Ziel ist die Wiedergabe einer außergewöhnlich breiten, räumlich präzisen und holografisch wirkenden Stereobühne über zwei Lautsprecher – insbesondere in Setups mit geringem Lautsprecherabstand.

Die ALSA-Loopback-Edition verwendet ein virtuelles ALSA-Audiogerät als Verbindung zwischen Audio-Player und DSP-Engine. Dadurch kann das Audiosignal systemintern abgegriffen, verlustfrei verarbeitet und anschließend nativ an den Hardware-DAC ausgegeben werden.

### Kerntechnologie: FFT-FIR & Overlap-Save

Die neueste FFT-FIR-Variante nutzt die extrem schnelle C-Bibliothek `libfftw3f` (Overlap-Save-Verfahren), um die CPU bei der Signalverarbeitung zu entlasten.

* **3-Wege-Frequenzweiche:** Das Signal wird im Frequenzbereich präzise aufgesplittet (Tief-Mittel-Hoch). Der RACE-Algorithmus wird ausschließlich auf den mittleren Frequenzbereich von empfohlenen ca. 200–400 Hz bis ca. 4500–5500 Hz beschränkt.
* **Phasenlinearer 3-Band Equalizer:** In der aktuellen Version wurden die klassischen IIR-Biquads restlos durch einen FFT-FIR-basierten Equalizer ersetzt. Die parametrischen EQ-Kurven werden rein als reelle Gain-Faktoren auf die komplexe FFT-Matrix multipliziert. Dies verhindert jegliche Phasenverschiebungen und garantiert eine perfekte Impulstreue (Transienten).
* **192 kHz / 32-Bit Float Pipeline:** Der DSP-Kern rechnet durchgehend mit 32-Bit Float bei 192 kHz. Integrierte XRUN-Recovery-Routinen und Clipping-Detektoren sorgen für höchste Stabilität im Langzeitbetrieb.

*Cave!* Bei schnellen Änderungen der Inputs über die Regler im Web-GUI kann es gelegentlich zu kleinen Störungen des Audio-Signals kommen, die jedoch innerhalb weniger Millisekunden wieder abklingen. Dies ist dem Algorithmus geschuldet, der den Signalfluss für die Neuberechnung der Arrays kurzzeitig unterbrechen muss. Die CPU-Last liegt bei korrekter Anwendung der Compiler-Flags (inkl. `-ffast-math -flto -fopenmp`) auf einem Raspberry Pi 5 trotz der extrem hohen Rechenintensität bei vertretbaren ca. 6 %.

Ist die Konfiguration für das jeweilige Setup einmal optimiert, können die zugehörigen Parameter in einer `oapw_config.txt` gespeichert werden. Sie werden beim nächsten Start automatisch eingelesen.

---

## 1. Systemanforderungen und Hardware

### Rechner

* **Raspberry Pi 5:** Empfohlene Referenz-Plattform für anspruchsvolle DSP-Verarbeitung mit nativen 192 kHz.
* **Raspberry Pi 4:** Alternative Plattform; die erreichbare maximale Abtastrate hängt von Konfiguration und Systemlast ab.
* **Raspberry Pi Compute Module:** Ebenfalls als Basis für dedizierte Audio-Systeme geeignet.

Der Raspberry Pi 3 wurde mit dieser Programmversion noch nicht abschließend auf seine Leistungsfähigkeit getestet.

### Audio-Hardware

Erforderlich ist ein hochwertiger, nativer 192-kHz-DAC (bevorzugt I²S-DAC-HAT), um das Phasenpotenzial des FIR-Algorithmus voll auszuschöpfen.

Mögliche Optionen sind:

* **Allo Katana** (Referenz-DAC)
* **Allo Boss DAC**
* **DAC HAT Pro**

*(Hinweis: Einfache USB-DACs, die auf 96 kHz limitiert sind (z. B. Sharkoon Gaming DAC), taugen nicht für diesen High-End-Anspruch und erzwingen ein verlustbehaftetes Downsampling über ALSA. Sie werden für die FIR-Edition ausdrücklich nicht empfohlen.)*

Die konkrete Konfiguration des ALSA-Gerätenamens (`playback_device`) in der `main.cpp` muss an die installierte Hardware angepasst werden.

### Stromversorgung

Für einen hochwertigen Audioaufbau ist eine sorgfältig ausgelegte, rauscharme Stromversorgung zwingend erforderlich (z. B. LPS - Linear Power Supply). Ein separates lineares Netzteil für den DAC reduziert Störeinflüsse massiv.

### WLAN und Bluetooth

WLAN und Bluetooth können bei einem dedizierten Audio-System deaktiviert werden, um Funkaktivität und unnötige Systemlast zu reduzieren.

---

## 2. Software-Voraussetzungen

Vorausgesetzt wird ein aktuelles 64-Bit-Debian-System, beispielsweise:

* Debian Trixie / Bookworm
* Raspberry Pi OS 64-Bit

Benötigte Komponenten:

* GCC/G++ mit C++17-Unterstützung
* GNU Make
* ALSA-Entwicklerbibliotheken (`libasound2-dev`)
* FFTW3-Bibliotheken (`libfftw3-dev`)
* Linux-Kernelmodul `snd-aloop`

### 2.1 Pakete installieren

```bash
sudo apt-get update
sudo apt-get install build-essential libasound2-dev libfftw3-dev

```

### 2.2 Repository herunterladen & Kompilieren

Das Repository (Zweig `alsa-edition`) klonen:

```bash
git clone -b alsa-edition https://github.com/druthibautgmail/OAPW_HiFi_FIR.git OAPW_HiFi_FIR
cd OAPW_HiFi_FIR

```

Das Programm kompilieren (nutzt alle verfügbaren Kerne für einen schnellen Build):

```bash
make -j4

```

Den Player für einen ersten Test manuell starten:

```bash
./OAPW_Player

```

### 2.3 WLAN und Bluetooth optional deaktivieren

Einträge in `/boot/firmware/config.txt` ergänzen:

```ini
dtoverlay=disable-wifi
dtoverlay=disable-bt

```

*(Bei älteren Installationen liegt die Datei unter `/boot/config.txt`.)*

### 2.4 ALSA-Loopback aktivieren

Das Kernelmodul `snd-aloop` stellt das virtuelle Audiokabel bereit.

Modul laden und Autostart aktivieren:

```bash
sudo modprobe snd-aloop
echo snd-aloop | sudo tee -a /etc/modules

```

Prüfen, ob das Modul geladen wurde:

```bash
lsmod | grep snd_aloop
aplay -l

```

---

## 3. Signalfluss – ALSA Loopback Pipeline (192 kHz / 32-Bit)

Der Signalweg verbindet den Audio-Player bitperfekt mit der OAPW-DSP-Engine.

```text
┌───────────────────────────────┐
│ Audio-Player (SOXR aktiv)     │
│ z. B. MPD oder myMPD          │
└───────────────┬───────────────┘
                │ 192 kHz
                ▼
┌───────────────────────────────┐
│ ALSA Loopback                 │
│ Playback: hw:Loopback,0,0     │
└───────────────┬───────────────┘
                │ Virtuelles Audiokabel
                ▼
┌───────────────────────────────┐
│ OAPW_HiFi_AL  (C++ Engine)    │
│ Capture: hw:Loopback,1,0      │
│                               │
│ • FFT Forward (libfftw3f)     │
│ • Phasenlinearer FIR-EQ       │
│ • RACE Crosstalk Elimination  │
│ • IFFT Backward               │
└───────────────┬───────────────┘
                │ 192 kHz / 32-Bit Float
                ▼
┌───────────────────────────────┐
│ Hardware-DAC                  │
│ z. B. Allo Katana (hw:Katana) │
└───────────────┬───────────────┘
                │
                ▼
┌───────────────────────────────┐
│ Verstärker und Lautsprecher   │
└───────────────────────────────┘

```

**Pegelstruktur (Gain Staging):** Um digitales Clipping bei der RACE-Matrix-Addition zu vermeiden, sollte der Audio-Player den Puffer nicht mit maximalem 0 dBFS Pegel beschreiben, oder es muss in der Web-GUI ein negativer *Makeup Gain* gesetzt werden.

---

## 4. Web-GUI und Steuerung

Die Engine startet standardmäßig einen lokalen Web-Server auf Port 8080.
Aufruf im Browser:

```text
http://<IP-Adresse-des-Raspberry-Pi>:8080

```

**Features der Web-GUI:**

* **Real-Time Analyzer:** Visualisierung des Echtzeit-Spektrums mit exaktem 192-kHz-Frequenz-Mapping und kalibrierter logarithmischer Dezibel-Skala.
* **Geometrie-Rechner:** Direkte Berechnung von Laufzeit (Delay) und Übersprechdämpfung (Attenuation) basierend auf dem realen Lautsprecher- und Hörerabstand.
* **DSP Parameter:** Steuerung der RACE-Matrix, Crossover-Frequenzen, des linearen FIR-Equalizers und des globalen Makeup-Gains.

---

## 5. Automatischer Start als Hintergrunddienst (systemd)

Damit `OAPW_Player` nach einem Neustart des Raspberry Pi automatisch hochfährt, wird er als `systemd`-Dienst eingerichtet.

### 5.1 Dienst-Datei erstellen

```bash
sudo nano /etc/systemd/system/oapw.service

```

Füge folgenden Block ein (Pfade und User anpassen!):

```ini
[Unit]
Description=OAPW HiFi ALSA Loopback Player
After=network.target sound.target mpd.service

[Service]
Type=simple
User=<Dein Username>
WorkingDirectory=/home/<Dein Username>/OAPW_HiFi_FIR
ExecStart=/home/<Dein Username>/OAPW_HiFi_FIR/OAPW_Player
StandardOutput=journal
StandardError=journal
Restart=always
RestartSec=5

[Install]
WantedBy=multi-user.target

```

### 5.2 Dienst aktivieren und verwalten

Dienst laden und starten:

```bash
sudo systemctl daemon-reload
sudo systemctl enable oapw.service
sudo systemctl start oapw.service

```

Logs zur Fehlerdiagnose (inkl. XRUNs und Clipping-Warnungen) abfragen:

```bash
sudo journalctl -u oapw.service -f

```

---

## 6. Hinweise zur Audioqualität und Performance

Die C++ Engine ist für kompromisslose Audioqualität bei 192 kHz ausgelegt. Wichtige Faktoren für den stabilen Betrieb:

* Pipewire und Wireplumber sollten systemweit deaktiviert/maskiert werden, um exklusiven ALSA-Zugriff zu garantieren.
* Der Zuspieler (z. B. MPD) muss das Quellmaterial zwingend vor der Übergabe an den Loopback resampeln (z. B. via SOXR), damit der DSP mit einem konstanten 192-kHz-Takt gefüttert wird.
* Thermisches Management: Der Raspberry Pi 5 muss ausreichend gekühlt werden, um Throttling zu vermeiden, da Einbrüche in der Taktrate sofort zu Puffer-Underruns (`-EPIPE`) auf ALSA-Ebene führen.

---

## 7. Danksagung

Mein ganz besonderer Dank gilt **Prof. Dr. Dieter Steinhilber**, dessen wertvolle Hinweise, technische Anregungen und tiefgreifendes Verständnis für High-End Audiophilie maßgeblich zur Weiterentwicklung und Reifung dieses Programms beigetragen haben.

Ein weiterer Dank gilt **Dr. Michael Graw**, meinem ehemaligen Ambiophonie-Mentor, der mir in meiner aktiven Zeit als CTO der Drägerwerk AG in Lübeck die Augen und vor allem die Ohren für die Ambiophonie geöffnet hat.

Ebenso gilt unser Dank **Ralph Glasgal** für sein Lebenswerk, seine wegweisenden Forschungen am Ambiophonics Institute und die Entwicklung der Ambiophonie, die dieses Softwareprojekt inspiriert und ermöglicht haben.

Dieses Programm wurde - insbesondere bei schnellen Änderungen und der Fehlersuche - mit Unterstützung durch Gemini und in kleinen Teilen auch ChatGPT entwickelt. Anregungen für das Design der Filter, des GUI und des Kern-Algorithmus stammen aus verschiedensten Quellen der Fachliteratur und der audiophilen Community, thanks and keep sharing!

---

*OAPW_HiFi_AL – Open Ambiophonics Processing Wrapper, High-End ALSA Loopback Edition.*

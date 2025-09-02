# My PC Hardware report

## Inxi package used to retrive info about PC hardware


Since the github repository seems to be **no longer actively maintained** [![Inxi-repo](https://git-scm.com/images/logos/downloads/Git-Logo-Black.png)](https://github.com/smxi/inxi), **it does it job well**. Einstein's *relativity theory* has stood the time until now, despite him not maintaining, if it works, don't fix it.

______________________________________________________________________

## 🚀 Features

Simply, it looks after variables such as:

* **System**
** *Subcategories such as:*
*** **Kernel: 6.16.4-arch1-1 arch: x86_64 bits: 64 compiler: gcc v: 15.2.1**
*** **Desktop: *Hyprland v: 0.50.1* Distro: *Arch Linux***
* **Machine**
* **Battey**
* **CPU**
* **GPU**
* **Audio**
* **Network**
* **Bluetooth**
* *And so on...*


## ⚙️ Installation

It may differ from distrubition to distribution (on arch linux, e.g.,):

```bash
sudo pacman -S inxi
```

Else, you find tarballs here:

[![Inxi-repo](https://upload.wikimedia.org/wikipedia/commons/9/9a/Codeberg_logo.svg)](https://github.com/smxi/inxi)

Just remember to:

```bash
tar -xvf file.tar.gz # <the downloaded tarball>
```

Running the command:

```bash
inxi -Fxz
```

### Example

Will produce an output like this:

```python
System:
  Kernel: 6.6.32-zen1-1-zen arch: x86_64 bits: 64 compiler: gcc v: 14.2.1
  Desktop: Hyprland v: 0.49.2 Distro: Arch Linux (rolling)
Machine:
  Type: Laptop System: Notebook (generic 14") v: 1.0 serial: [redacted]
  Mobo: Generic model: MB14 v: 1.0 serial: [redacted]
    UEFI: American Megatrends v: 1.14 date: 2022-06-21
Battery:
  ID-1: BAT0 charge: 41.3 Wh (82.6%) condition: 49.8/52.0 Wh (95.7%)
    volts: 12.1 min: 11.4 model: Standard Battery status: discharging
CPU:
  Info: 6-core model: AMD Ryzen 5 5500U with Radeon Graphics
    bits: 64 type: MT MCP arch: Zen 2 rev: 1 cache: L1: 384 KiB
    L2: 3 MiB L3: 8 MiB
  Speed (MHz): avg: 2100 min/max: 400/4050 boost: enabled cores: 1: 2100
    2: 2100 3: 2100 4: 2100 5: 2100 6: 2100 7: 2100 8: 2100 bogomips: 35900
  Flags: avx avx2 ht lm nx pae sse sse2 sse3 sse4_1 sse4_2 sse4a ssse3 svm
Graphics:
  Device-1: Advanced Micro Devices [AMD/ATI] Lucienne vendor: Generic
    driver: amdgpu v: kernel arch: GCN-5 bus-ID: 03:00.0 temp: 46.0 C
  Device-2: USB HD UVC WebCam driver: uvcvideo type: USB bus-ID: 3-2:2
  Display: wayland server: X.org v: 1.21.1.18 with: Xwayland v: 24.1.8
    compositor: Hyprland v: 0.49.2 driver: X: loaded: amdgpu
    unloaded: modesetting dri: radeonsi gpu: amdgpu resolution: 1920x1080~60Hz
  API: Vulkan v: 1.3.280 drivers: radv devices: 1
  API: EGL Message: EGL data requires eglinfo. Check --recommends.
  Info: Tools: api: vulkaninfo wl: nwg-displays x11: xprop,xrandr
Audio:
  Device-1: AMD Raven/Renoir/Vangogh HDMI/DP Audio driver: snd_hda_intel
    v: kernel bus-ID: 03:00.1
  Device-2: AMD Family 17h/19h HD Audio driver: snd_hda_intel v: kernel
    bus-ID: 03:00.6
  API: ALSA v: k6.6.32-zen1-1-zen status: kernel-api
  Server-1: JACK v: 1.9.22 status: off
  Server-2: PipeWire v: 1.2.1 status: active
Network:
  Device-1: Intel Wi-Fi 6 AX200 driver: iwlwifi v: kernel bus-ID: 02:00.0
  IF: wlan0 state: up mac: [hidden]
Bluetooth:
  Device-1: Intel AX200 Bluetooth driver: btusb v: 0.8 type: USB bus-ID: 3-1:2
  Report: btmgmt ID: hci0 rfk-id: 2 state: up address: [hidden]
Drives:
  Local Storage: total: 476.94 GiB used: 182.11 GiB (38.2%)
  ID-1: /dev/nvme0n1 vendor: [generic] model: NVMe SSD size: 476.94 GiB
    temp: 42.0 C
Partition:
  ID-1: / size: 96.1 GiB used: 61.7 GiB (64.2%) fs: ext4 dev: /dev/dm-0
    mapped: vg0-root
  ID-2: /boot size: 1.0 GiB used: 286.0 MiB (27.5%) fs: ext4
    dev: /dev/nvme0n1p2
  ID-3: /home size: 128.3 GiB used: 70.2 GiB (54.7%) fs: ext4
    dev: /dev/dm-1 mapped: vg0-home
Swap:
  ID-1: swap-1 type: zram size: 8 GiB used: 0 KiB (0.0%) dev: /dev/zram0
  ID-2: swap-2 type: partition size: 8 GiB used: 0 KiB (0.0%)
    dev: /dev/nvme0n1p3
Sensors:
  System Temperatures: cpu: 46.8 C gpu: amdgpu temp: 45.0 C
  Fan Speeds (rpm): cpu: 2200
Info:
  Memory: total: 16 GiB available: 12.3 GiB used: 3.1 GiB (25.2%)
  Processes: 245 Uptime: 1h 12m Init: systemd
  Packages: 2860 Compilers: clang: 18.1.8 gcc: 14.2.1 Shell: Zsh v: 5.9
    inxi: 3.3.36
```


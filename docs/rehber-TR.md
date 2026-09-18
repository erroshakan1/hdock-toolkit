# HDOCK rehberi — mantığıyla birlikte

Bu alana yeni girmiş biri için. Her adımda **ne yaptığımız** değil, **neden ona
baktığımız** anlatılıyor.

---

## 1. Docking ne yapar, ne yapmaz

İki molekülü alır, birini sabit tutup (reseptör) diğerini (ligand) bütün göreli
konum ve açılarda dener, her denemeyi puanlar. FFT kullandığı için milyonlarca
konumu dakikalar içinde tarar.

**Yapmadığı şey:** moleküllerin şekil değiştirmesine izin vermez. İkisi de rijit
gövde. Çözücü (su) açıkça yok, örtük. Yani docking sana "bunlar bağlanır"
demiyor; "bağlansalardı geometrik olarak en uygun duruş bu olurdu" diyor.

Sonuç bir **başlangıç yapısıdır**, bir sonuç değil. Asıl testi MD yapar.

---

## 2. Reseptör mü ligand mı

Büyük olan reseptör, küçük olan ligand. Sebebi algoritmik: büyük olanın etrafında
ızgara kurulup küçük olan o ızgarada gezdiriliyor.

---

## 3. Girdi hazırlığı neden gerekli

```bash
python3 scripts/prepare_pdb.py giris.pdb cikis.pdb [zincir]
```

Üç şey yapar:

- **Hidrojenleri atar.** Docking ağır atom geometrisiyle çalışır; hidrojenler
  gereksiz hesap yükü ve bazı programlarda ayrıştırma hatası sebebi.
- **Zincir kimliği yazar.** `gmx`'in ürettiği PDB'lerde 22. sütun boş kalır.
  Zincir kimliği olmadan ligandı reseptörden ayıramazsın.
- **Element sütununu doldurur** (77–78. sütun). Bazı araçlar buradan atom tipini
  okur; boşsa yanlış yorumlar.

> **Tuzak.** Ligandın zincir harfi reseptörünkilerle çakışmamalı. Çakışırsa model
> dosyasında ikisi tek zincirde birleşir. Skorlar doğru kalır ama **kümeleme
> analizi sahte sonuç verir**: ligandın kütle merkezi diye hesapladığın şeyin
> içine hiç kıpırdamayan reseptör atomları girer, bütün pozlar birbirine yakın
> görünür. Bu projede bir kez yaşandı — S6'da "9/9 poz toplandı" çıktı, gerçeği
> 0–5/9'du. `run_docking.sh` artık boş harfi otomatik seçiyor.

---

## 4. Skorlar

### Docking skoru
Birimsiz. **Daha negatif = daha iyi.** Enerji değil, istatistiksel bilgi-tabanlı
bir potansiyel. kJ/mol'e çevrilemez, **asla ΔG olarak rapor edilmez.** Sadece
sıralama için.

### Güven skoru
Docking skorundan hesaplanan kalibre olasılık:

```
güven = 1 / (1 + exp(0,02 × (skor + 150)))
```

| Güven | Yorum |
|---|---|
| > 0,7 | bağlanma çok muhtemel |
| 0,5–0,7 | olası |
| < 0,5 | zayıf |

Kalibrasyon protein–protein kıyaslama setlerinde yapıldı. Protein–DNA için yön
gösterir, kanıt değildir.

### Ligand RMSD — yanıltıcı sütun
Pozun **girdi olarak verdiğin** ligand koordinatlarından uzaklığı. Girdiyi
rastgele yerleştirdiysen anlamsızdır; model 1'in RMSD'si bile sıfır çıkmaz.
Yalnızca girdi *bilinen bir kompleks* ise anlamlıdır — o zaman düşük RMSD
"docking doğru cevabı buldu" demektir (bkz. `examples/1cgi.md`).

---

## 5. Asıl bakılan şey: pozlar tekrar ediyor mu

Tek bir iyi skor gürültü olabilir. **Farklı ligand konformerlerinden gelen
pozların aynı bölgeyi bulması** olamaz.

```bash
python3 scripts/select_pose.py docking/isim
```

Her poz için reseptörün hangi kalıntılarına 4,5 Å'den yakın olduğu çıkarılır —
buna **temas parmak izi** denir. İki pozun parmak izleri %50'den fazla örtüşüyorsa
aynı kümeye girerler (Jaccard benzerliği). En kalabalık kümenin en iyi skorlusu
başlangıç yapısı olur.

---

## 6. Kalıntı düzeyinde harita

```bash
python3 scripts/contact_map.py docking/isim --top 15
```

Her kalıntı için "kaç pozda temas edildi" sayar. **Poz başına bir kez sayar**,
atom çifti başına değil — bu ayrımı kaçırırsan yüzdeler 900'e çıkar.

Çok pozda görünen kalıntı gerçek sıcak nokta; bir pozda görünen gürültü.

---

## 7. α-sarmal yüzü — asıl anlatmak istediğim kısım

### Kalıntı ne demek
Protein bir zincir, halkaları amino asitler. Her amino asit bir **kalıntı**
(residue) ve zincir boyunca numaralanır: 1, 2, 3… "TFAM'ın 60. kalıntısı"
demek, zincirin başından 60. amino asit demek.

### Neden numaraların *aralığına* bakıyoruz

Bir α-sarmalda tam **3,6 kalıntı bir tur** atar. Buradan:

```
360° ÷ 3,6 kalıntı = kalıntı başına 100° dönüş
```

Zincirde ilerlerken her kalıntı bir öncekinden 100° dönmüş olur. Sonucu:

| Kaç adım ileri | Toplam dönüş | Nereye bakar |
|---|---|---|
| +1 | 100° | tamamen başka yön |
| +2 | 200° | karşı taraf |
| **+3** | 300° | ≈ aynı yön |
| **+4** | 400° → 40° | ≈ aynı yön |
| +7 | 700° → 340° | ≈ aynı yön |

Yani **i, i+3, i+4, i+7** kalıntıları sarmalın **aynı yüzünde** toplanır.
"Amfipatik sarmal", "heptad tekrarı" gibi kavramların tamamı bu geometriden
çıkar.

### Somut örnek

Aβ42'nin TFAM üzerinde temas ettiği kalıntılar:

```
60 · 63 · 64 · 67 · 68 · 71 · 72
   +3   +1   +3   +1   +3   +1
```

i, i+3, i+4 örüntüsü. Demek ki hepsi tek bir yüzde.

### Doğrulaması — sarmal tekerleği

```bash
python3 scripts/helix_wheel.py model.pdb --chain A --range 56-71 \
        --contacts 60,63,64,67,68,71 --offset 43
```

Sarmalın eksenini uydurur, her CA'yı eksene dik düzleme yansıtır, açısını verir:

| kalıntı | açı | temas |
|---|---|---|
| 60 | 321° | ✓ |
| 61 | 87° | |
| 62 | 197° | |
| 63 | 264° | ✓ |
| 64 | 328° | ✓ |
| 67 | 293° | ✓ |
| 68 | 41° | ✓ |
| 71 | 10° | ✓ |

Temas edenler dar bir yayda (dairesel yayılım 47°), temas etmeyenler sarmalın
öbür taraflarında. **Kalıntı başına ölçülen dönüş 89°** — teorik α-sarmal
değeri 100°, tutuyor.

### Bu ne eliyor

- **Rastgele yüzey teması olsaydı** açılar 0–360° arasına yayılırdı → docking
  ligandı gelişigüzel bir yere kondurmuş olurdu.
- **Ardışık kalıntılar (60, 61, 62, 63) çıksaydı** bu geometrik olarak imkânsızdı
  — sarmalın her yönüne aynı anda dokunamazsın. Öyle bir liste halkaya veya
  zincir ucuna işaret ederdi.

Kalan tek açıklama: tanımlı, yönlü bir yüzey yaması. "Bir cep var" diyebiliyoruz.

### Tuzak — aralığı gözle seçme

Sarmal sınırını DSSP'ye sor, tahmin etme. Bu projede 58–76 aralığı denendi, ama
DSSP'ye göre 72–76 halkaydı. Halkayı içeren aralıkta eksen uydurması bozuldu ve
kalıntı başına dönme 73° çıktı. Gerçek sarmal sınırıyla (56–71) 89° oldu.

```python
from MDAnalysis.analysis.dssp import DSSP
print("".join(DSSP(u.select_atoms("protein")).run().results.dssp[0]))
# H = alfa-sarmal, E = beta-tabaka, - = halka
```

---

## 8. Sınırlar — ne iddia edemezsin

- Rijit gövde: ligand da reseptör de esnemiyor. Düzensiz peptitler için özellikle
  kısıtlayıcı; biz bu yüzden MD'den gelen **üç farklı konformer** kullandık.
- Çözücü örtük: spesifik su köprüleri ve iyonlar yok.
- Skor mutlak afinite değil.
- Güven skoru protein–protein kalibrasyonu.

Docking'in çıktısı **MD'nin başlangıç yapısıdır**. Makalede "docking gösterdi ki"
diye cümle kurma; "docking başlangıç yapısını üretti, MD şunu gösterdi" de.

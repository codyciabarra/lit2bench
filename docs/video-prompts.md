# Higgsfield prompt set — Lit2Bench site footage

Generation prompts for the atmospheric video and stills on the Lit2Bench site.
Everything here is **decorative b-roll**: light, glass, fluid, lab texture.

**Nothing in this file generates product footage, and it shouldn't.** A model
asked for "the Lit2Bench interface" will invent a UI that isn't ours and, far
worse, invent sashimi plots, junction counts and ΔΔCt values. Fabricated results
on the landing page of a cryptic-splicing tool contradicts the exact claim the
site is built on. Product demos are screen recordings of the real app.

That is also why every prompt below carries `no text, no UI, no charts, no
graphs` — left unsaid, these models garnish scientific scenes with garbled
labels and plausible-looking plots.

---

## Brand lock

Append to **every** prompt so the clips cut together as one piece. Same palette
as `site/styles.css` and the app's `--l2b-*` tokens.

```
Shot on ARRI Alexa, 35mm anamorphic, shallow depth of field, 24fps.
Palette strictly violet #8f7dfa, electric blue #4f8cff and teal #2fd9c4
against near-black #080a13. Volumetric haze, soft rim light, single key
source, deep falloff into pure black. High-end pharmaceutical commercial
aesthetic, clinical and restrained, not neon and not cyberpunk.
No text, no logos, no UI, no charts, no graphs, no on-screen graphics.
```

**Negative prompt** (paste into the negative field every time):

```
text, letters, numbers, watermark, logo, subtitles, user interface, charts,
graphs, plots, diagrams, scientific figures, distorted hands, extra fingers,
malformed glassware, oversaturated, neon, cyberpunk, cartoon, 3D render look,
plastic skin, dirty lens, heavy lens flare
```

Lock a **seed** on the first clip you like and reuse it across the set. Palette
drift between sections is the single most obvious tell that footage was
generated piecemeal.

---

## 1. Hero ambient loop — the important one

Replaces or sits behind the CSS blobs in `.hero`. This is the clip people
actually see, so generate 6–8 variations and keep one.

- **Model:** DoP / text-to-video · **Preset:** `Dolly In` at minimum speed (or `Static`)
- **16:9 · 5s · seamless loop**

```
Extreme macro of a translucent glass double-helix ribbon suspended and slowly
rotating in perfectly still black liquid. Light refracts through the glass into
violet, electric blue and teal. Fine luminous particles drift slowly past the
lens. Nothing else in frame, absolute darkness beyond the subject.
```

Motion must be **almost imperceptible** — it plays under a headline. Anything
faster fights the type instead of supporting it.

---

## 2. "Detect" — signal separating from noise

For the Workflow section, step 01. The one clip that carries real meaning: a
single strand resolving out of an indistinguishable mass.

- **Preset:** `Focus Change` (rack focus) · **16:9 · 5s**

```
Extreme macro through a shallow layer of clear fluid. Hundreds of fine parallel
filaments of light drift in soft focus. One filament brightens to teal and
separates from the mass as focus racks slowly from the blurred crowd onto that
single sharp strand.
```

---

## 3. "Design" — the pipette

Step 02.

- **Preset:** `Dolly In` or `Robo Arm` · **16:9 · 4s**

```
Extreme macro of a precision pipette tip releasing one droplet toward a well
plate. The droplet catches violet and teal light in mid-fall, surface tension
deforming it. Black lab bench, hard rim light along the pipette edge, deep
shadow everywhere else.
```

---

## 4. "Validate" — the plate

Step 03. Keep the well glow **graded and even** — a pattern of bright and dark
wells reads as data, which is exactly what this must not do.

- **Preset:** `Crane Down` · **16:9 · 4s**

```
Top-down macro of a clear 96-well plate on black glass, wells holding an even
gradient of teal into violet, faint condensation on the seal. Slow crane down
toward the surface. Soft even overhead light.
```

---

## 5. "Record" — the notebook

Step 04. The one warm, human shot in the set; it stops the sequence feeling
like a screensaver.

- **Preset:** `Static` or gentle `Handheld` · **16:9 · 4s**

```
Overhead macro of a hand in a blue nitrile glove closing a worn hardcover lab
notebook on a dark bench, a rack of capped tubes just out of focus behind it.
Single warm key light from the left, shallow depth of field, blank unmarked
pages.
```

`blank unmarked pages` is doing real work — without it you get pages of
hallucinated handwriting.

---

## 6. Neuron — for the Gitler Lab credits section

Grounds the site in the actual biology the toolkit targets.

- **Preset:** `360 Orbit` (slow) or `Push In` · **16:9 · 6s**

```
Extreme macro of a single neuron rendered as translucent blown glass, dendrites
branching away into darkness, faint teal pulses travelling along the axon. Slow
orbit around the cell body. Suspended in black, dust motes catching the light.
```

---

## 7. Open-graph still + video posters

Every `<video>` needs a `poster`, or the section flashes empty on slow
connections. Generate these as **stills**, not frames grabbed from the clips —
grabs are motion-blurred.

- **Model:** Soul · **16:9 · 1920×1080**

```
Extreme macro still life: a translucent glass double helix standing on black
glass, refracting violet, blue and teal light, volumetric haze behind it.
Centred, generous negative space in the upper third. Editorial product
photography, ARRI Alexa, 35mm anamorphic.
```

Negative space up top is deliberate — that's where the title sits in a link
preview card.

## 8. Social cuts

Re-render clips 1, 2 and 6 at **9:16** rather than cropping. Anamorphic framing
crops badly — the subject ends up half out of frame.

---

## Putting them on the site

Constraints the current build imposes:

- **Budget ~2.5 MB per loop.** 1920×1080, 24fps, 5s, H.264 CRF 24 lands there.
  Ship WebM/VP9 too; Safari takes the MP4, everything else the smaller WebM.
- **Everything stays local.** No CDN — the site makes zero external requests and
  that should hold.
- **`preload="none"` on anything below the fold**, or six clips download on load
  and the page feels slower than it did with no video at all.

```html
<video class="hero-video" autoplay muted loop playsinline
       poster="assets/hero-poster.jpg" preload="none" aria-hidden="true">
  <source src="assets/hero.webm" type="video/webm">
  <source src="assets/hero.mp4"  type="video/mp4">
</video>
```

- **`prefers-reduced-motion` must show the poster only.** The site already
  honours that setting everywhere else; an autoplaying video would be the one
  thing that ignores it.

```css
@media (prefers-reduced-motion: reduce) {
  .hero-video { display: none; }
  .hero-media { background: url("assets/hero-poster.jpg") center/cover; }
}
```

- `aria-hidden="true"` and no audio track: this is decoration, and a screen
  reader has nothing to say about it.

Encoding, once you have the renders:

```bash
ffmpeg -i hero-raw.mp4 -an -c:v libx264 -crf 24 -preset slow \
       -pix_fmt yuv420p -movflags +faststart assets/hero.mp4
ffmpeg -i hero-raw.mp4 -an -c:v libvpx-vp9 -crf 34 -b:v 0 assets/hero.webm
ffmpeg -i hero-raw.mp4 -vframes 1 -q:v 3 assets/hero-poster.jpg
```

## The real product demo

Screen-record the app. The sequence worth capturing, in one take:

1. Cryptic Engine — load a BAM pair, run detection, sashimi plot appears
2. Double-click to zoom the plot
3. Click a candidate exon → "Design primers for this exon"
4. Primer designer opens pre-filled; the schematic renders

That hand-off between tools is the actual product story, and no generated clip
can tell it. `⇧⌘5` records a window; keep it under 30s, no cursor trails, and
crop to the window with a couple of pixels of margin.

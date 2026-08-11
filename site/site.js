/* site.js -- interaction for the Lit2Bench site.
 *
 * No framework and no bundler: this is served straight off GitHub Pages next to
 * index.html. Everything here is progressive -- the page reads and works with
 * JavaScript switched off, and each block below only enhances what's already
 * rendered in the HTML.
 */
(function () {
  'use strict';

  /* ------------------------------------------------------------- release
   * The one place download URLs are defined. Point this at whatever the
   * current release is; everything on the page follows from it.
   *
   * `latest/download/<asset>` is a GitHub redirect that always resolves to the
   * newest release's asset, so this keeps working across releases without an
   * edit here -- but the asset filename carries the version, so `version` and
   * `asset` do have to move together.
   */
  var RELEASE = {
    version: '0.1.0',
    size: '~1.4 MB',
    repo: 'https://github.com/codyciabarra/lit2bench',
    asset: 'Lit2Bench-0.1.0.dmg'
  };
  RELEASE.url = RELEASE.repo + '/releases/latest/download/' + RELEASE.asset;

  var reduced = matchMedia('(prefers-reduced-motion: reduce)').matches;
  var $  = function (s, r) { return (r || document).querySelector(s); };
  var $$ = function (s, r) { return Array.prototype.slice.call((r || document).querySelectorAll(s)); };

  /* --------------------------------------------------------------- download */
  var dlBtn = $('#dl-btn');
  if (dlBtn) {
    dlBtn.href = RELEASE.url;
    $('#dl-version').textContent = RELEASE.version;
    $('#dl-size').textContent = RELEASE.size;

    // Visitors on Windows or Linux shouldn't be handed a .dmg without warning:
    // the app runs fine from source on their machine, only the installer is
    // macOS-only so far.
    var ua = navigator.userAgent;
    if (!/Mac/i.test(ua) || /iPhone|iPad/i.test(ua)) {
      $('#dl-label').textContent = 'Download the macOS .dmg';
      $('#dl-sub').textContent =
        'macOS only for now — on Windows or Linux, run from source (see right).';
    }
  }
  var year = $('#year');
  if (year) year.textContent = new Date().getFullYear();

  /* ------------------------------------------------------------------- nav */
  var nav = $('#nav');
  var onScroll = function () { nav.classList.toggle('scrolled', window.scrollY > 8); };
  onScroll();
  addEventListener('scroll', onScroll, { passive: true });

  /* ----------------------------------------------------------------- theme */
  var themeBtn = $('#theme');
  if (themeBtn) {
    themeBtn.addEventListener('click', function () {
      var next = document.documentElement.dataset.theme === 'light' ? 'dark' : 'light';
      document.documentElement.dataset.theme = next;
      try { localStorage.setItem('l2b-theme', next); } catch (e) {}
    });
  }

  /* ---------------------------------------------------------------- reveal */
  var revealables = $$('.reveal');
  if (reduced || !('IntersectionObserver' in window)) {
    revealables.forEach(function (el) { el.classList.add('in'); });
  } else {
    var io = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        if (!e.isIntersecting) return;
        // Stagger siblings so a row of cards arrives as a wave, not a slab.
        var peers = $$('.reveal', e.target.parentNode);
        e.target.style.transitionDelay = (Math.min(peers.indexOf(e.target), 6) * 60) + 'ms';
        e.target.classList.add('in');
        io.unobserve(e.target);
      });
    }, { rootMargin: '0px 0px -12% 0px', threshold: 0.08 });
    revealables.forEach(function (el) { io.observe(el); });

    // Belt and braces: some renderers (headless Chrome without a compositor,
    // and a few in-app webviews) never deliver intersection callbacks. Reveal
    // everything unconditionally after a moment so a missed callback can't
    // leave a section permanently invisible.
    setTimeout(function () {
      revealables.forEach(function (el) { el.classList.add('in'); });
    }, 2500);
  }

  /* -------------------------------------------------------------- counters */
  // Numbers count up once, the first time they scroll into view. The element's
  // text is already correct in the HTML, so this is decoration only.
  $$('[data-count]').forEach(function (el) {
    var target = parseInt(el.dataset.count, 10);
    if (reduced || isNaN(target) || target === 0 || !('IntersectionObserver' in window)) return;
    var obs = new IntersectionObserver(function (entries) {
      if (!entries[0].isIntersecting) return;
      obs.disconnect();
      var t0 = performance.now(), dur = 900;
      (function tick(now) {
        var p = Math.min((now - t0) / dur, 1);
        el.textContent = Math.round(target * (1 - Math.pow(1 - p, 3)));
        if (p < 1) requestAnimationFrame(tick);
      })(t0);
    }, { threshold: 0.6 });
    obs.observe(el);
  });

  /* ---------------------------------------------------------- tool filters */
  var grid = $('#tool-grid');
  if (grid) {
    var tools = $$('.tool', grid);
    var chips = $$('.chip');
    var search = $('#tool-search');
    var empty = $('#no-results');
    var group = 'all';

    var apply = function () {
      var q = (search.value || '').trim().toLowerCase();
      var shown = 0;
      tools.forEach(function (t) {
        var okGroup = group === 'all' || t.dataset.group === group;
        var okText = !q || t.textContent.toLowerCase().indexOf(q) !== -1;
        var show = okGroup && okText;
        t.classList.toggle('hidden', !show);
        if (show) shown++;
      });
      empty.hidden = shown > 0;
    };

    chips.forEach(function (chip) {
      chip.addEventListener('click', function () {
        group = chip.dataset.filter;
        chips.forEach(function (c) { c.setAttribute('aria-pressed', String(c === chip)); });
        apply();
      });
    });
    search.addEventListener('input', apply);

    /* Pointer spotlight: feed the cursor position to the card's ::before
       gradient. Coordinates are per-card and only updated while hovering, so
       this costs nothing when the pointer is elsewhere. */
    if (!reduced && matchMedia('(hover: hover)').matches) {
      tools.forEach(function (t) {
        t.addEventListener('pointermove', function (e) {
          var r = t.getBoundingClientRect();
          t.style.setProperty('--mx', (e.clientX - r.left) + 'px');
          t.style.setProperty('--my', (e.clientY - r.top) + 'px');
        });
      });
    }
  }

  /* ------------------------------------------------------------------ copy */
  $$('.copy').forEach(function (btn) {
    btn.addEventListener('click', function () {
      var code = $('code', btn.parentNode);
      if (!code || !navigator.clipboard) return;
      navigator.clipboard.writeText(code.textContent).then(function () {
        var was = btn.textContent;
        btn.textContent = 'Copied';
        btn.classList.add('ok');
        setTimeout(function () { btn.textContent = was; btn.classList.remove('ok'); }, 1600);
      });
    });
  });

  /* ----------------------------------------------------------- b-roll clips */
  // The decorative clips carry no `autoplay` attribute: they start when they
  // scroll into view and pause again on the way out, so five of them never sit
  // off-screen decoding at once. With reduced motion, or without an
  // IntersectionObserver, nothing is ever started and each <video> shows its
  // poster instead -- which is exactly the intended fallback, so there is no
  // separate no-JS path to maintain.
  var brolls = $$('[data-broll]');
  if (brolls.length && !reduced && 'IntersectionObserver' in window) {
    var bio = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        var v = e.target;
        if (!e.isIntersecting) { if (!v.paused) v.pause(); return; }
        // preload="none" means nothing is buffered until we ask, so load()
        // first -- otherwise the first play() shows a visible stall.
        if (!v.dataset.loaded) {
          // Not every clip in the set has been shot yet. A missing file would
          // otherwise leave an empty bordered box on the page, so a clip that
          // cannot load takes its container down with it and the layout closes
          // up as though it had never been there.
          //
          // Two mechanisms, because neither is sufficient alone. A failing
          // <source> fires `error` at the <source> itself and `error` does not
          // bubble -- so the listener has to run in the *capture* phase to see
          // it from up here. And when the browser exhausts every <source> it
          // sets NETWORK_NO_SOURCE without firing anything at the <video> at
          // all, which only a check after the fact catches.
          var hide = function () {
            var box = v.closest('.broll');
            if (box) box.hidden = true;
          };
          v.addEventListener('error', hide, { capture: true, once: true });
          v.load();
          setTimeout(function () {
            if (v.networkState === 3 /* NETWORK_NO_SOURCE */) hide();
          }, 1200);
          v.dataset.loaded = '1';
        }
        var playing = v.play();
        // Autoplay can still be refused (low power mode, a strict setting).
        // That leaves the poster up, which is a perfectly good outcome.
        if (playing && playing.catch) playing.catch(function () {});
      });
    }, { threshold: 0.25 });
    brolls.forEach(function (v) { bio.observe(v); });
  }

  // The hero clip is the one that autoplays from the markup, so reduced motion
  // has to stop it here too: `display:none` hides a video without reliably
  // stopping it decoding.
  var heroVideo = $('.hero-video');
  if (heroVideo && reduced) {
    heroVideo.removeAttribute('autoplay');
    heroVideo.pause();
  }

  /* ------------------------------------------------------- screenshot tilt */
  // A few degrees of parallax on the hero shot as the pointer crosses the
  // viewport. Skipped entirely on touch and under reduced motion.
  var tilt = $('.shot-tilt');
  if (tilt && !reduced && matchMedia('(hover: hover)').matches) {
    addEventListener('pointermove', function (e) {
      var x = (e.clientX / innerWidth - 0.5) * 2;
      var y = (e.clientY / innerHeight - 0.5) * 2;
      tilt.style.setProperty('--tilt-y', (x * 2.2).toFixed(2) + 'deg');
      tilt.style.setProperty('--tilt-x', (-y * 1.4).toFixed(2) + 'deg');
    }, { passive: true });
  }
})();

---
layout: gallery
title: "Gallery"
permalink: /gallery/
excerpt: "A visual field notebook of travel, academic journeys, architecture, landscape, and quiet observation."
author_profile: false
hide_title: true
---

<div class="artist-gallery">
  <section class="artist-gallery__hero" aria-labelledby="gallery-title">
    <div class="artist-gallery__hero-copy">
      <p class="artist-gallery__eyebrow">A VISUAL FIELD NOTEBOOK · 2022—2026</p>
      <h1 id="gallery-title">Light, distance,<br>and the spaces between.</h1>
      <p class="artist-gallery__lede">A personal archive gathered between research journeys: mountains emerging from cloud, cities after dark, rooms made quiet by art, and brief encounters with unfamiliar places.</p>
      <a class="artist-gallery__scroll" href="#gallery-archive">Explore the albums <span aria-hidden="true">↓</span></a>
    </div>

    <div class="artist-gallery__poster" aria-label="Selected photographs">
      <figure class="artist-gallery__poster-main">
        <img src="{{ '/images/gallery/2025-06-07_Nepal_Himalayan_Valley_0804.jpg' | relative_url }}" alt="Himalayan ridges and valleys in Nepal" width="1800" height="1350" fetchpriority="high">
        <figcaption>Nepal · 2025</figcaption>
      </figure>
      <figure class="artist-gallery__poster-top">
        <img src="{{ '/images/gallery/2024-07-14_London_UK_Contemporary_Gallery_1513.jpg' | relative_url }}" alt="A quiet contemporary gallery in London" width="1800" height="1350" fetchpriority="high">
        <figcaption>London · 2024</figcaption>
      </figure>
      <figure class="artist-gallery__poster-bottom">
        <img src="{{ '/images/gallery/2025-12-22_Nevada_USA_Red_Rock_Portrait_1353.jpg' | relative_url }}" alt="Portrait among the red rocks of Nevada" width="1350" height="1800" fetchpriority="high">
        <figcaption>Nevada · 2025</figcaption>
      </figure>
      <span class="artist-gallery__poster-mark" aria-hidden="true">FIELD<br>NOTES</span>
    </div>
  </section>

  <section class="artist-gallery__statement" aria-label="Artist statement">
    <span class="artist-gallery__statement-number">01</span>
    <p>These are not records of arrival. They are studies of atmosphere—how a place holds light, how scale changes attention, and how travel turns ordinary moments into composition.</p>
  </section>

  <section class="artist-gallery__archive" id="gallery-archive" aria-labelledby="archive-title">
    <header class="artist-gallery__archive-head">
      <div>
        <p class="artist-gallery__eyebrow">THE ARCHIVE</p>
        <h2 id="archive-title">Collected observations</h2>
      </div>
      <p class="artist-gallery__archive-note">{{ site.data.gallery_albums | size }} albums · {{ site.data.gallery | size }} photographs</p>
    </header>

    <div class="artist-gallery__albums">
      {% for album in site.data.gallery_albums %}
        {% assign album_url = '/gallery/' | append: album.slug | append: '/' %}
        {% assign cover_url = '/images/gallery/covers/' | append: album.slug | append: '.webp' %}
        {% assign cover = site.data.gallery | where: "id", album.cover_id | first %}
        <a class="artist-gallery__album" href="{{ album_url | relative_url }}" aria-labelledby="album-{{ album.slug }}">
          <div class="artist-gallery__album-cover">
            <img src="{{ cover_url | relative_url }}" alt="{{ cover.alt | escape }}" width="900" height="600" loading="lazy" decoding="async">
          </div>
          <div class="artist-gallery__album-copy">
            <p class="artist-gallery__album-count">{{ album.photo_ids | size }} photographs</p>
            <h3 id="album-{{ album.slug }}">{{ album.title | escape }}</h3>
            <p class="artist-gallery__album-description">{{ album.description | escape }}</p>
            <span class="artist-gallery__album-arrow" aria-hidden="true">↗</span>
          </div>
        </a>
      {% endfor %}
    </div>
  </section>

  <section class="artist-gallery__closing" aria-label="Closing note">
    <p class="artist-gallery__eyebrow">A CONTINUING PRACTICE</p>
    <p class="artist-gallery__closing-title">The archive remains open.</p>
    <p>New places will enter slowly—one frame, one encounter, one remembered light at a time.</p>
    <a href="{{ '/amateur/' | relative_url }}">Continue to film and visual storytelling <span aria-hidden="true">↗</span></a>
  </section>
</div>

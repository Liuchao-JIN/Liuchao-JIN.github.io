# 个人主页优化方案

> 针对 https://liuchao-jin.github.io/ （AcademicPages / Minimal Mistakes based Jekyll 站点）
> 目标：**用户友好、美观、内容丰富**
> 整理日期：2026-07-17

---

## 目录

- [一、总体诊断](#一总体诊断)
- [二、优先级总览（先做哪些）](#二优先级总览先做哪些)
- [三、首页（About）改造](#三首页about改造)
- [四、新增 News / 动态板块](#四新增-news--动态板块)
- [五、Publications 页面升级](#五publications-页面升级)
- [六、CV 页面结构化](#六cv-页面结构化)
- [七、视觉与美观（主题 / 配色 / 排版）](#七视觉与美观主题--配色--排版)
- [八、导航与信息架构](#八导航与信息架构)
- [九、SEO 与可被搜索性](#九seo-与可被搜索性)
- [十、性能与可访问性](#十性能与可访问性)
- [十一、技术债与链接规范](#十一技术债与链接规范)
- [十二、加分项（锦上添花）](#十二加分项锦上添花)
- [十三、落地清单（Checklist）](#十三落地清单checklist)

---

## 一、总体诊断

**优点**
- 内容非常扎实：13 篇论文、专利、大量教学记录、丰富的获奖与评审服务、甚至业余影视/阅读作品。
- 已配置好 ORCID、Google Scholar、ResearchGate 等学术身份链接。
- 使用成熟的学术主页模板，结构清晰、可维护。

**主要问题（可优化点）**
| 问题 | 现状 | 影响 |
|---|---|---|
| 首页信息量低 | About 只有 1 段文字 + 1 张封面图 | 第一印象弱，无法快速展示成果 |
| 缺少 News/动态 | 无 | 访客无法感知"活跃度"，学术主页几乎必备 |
| 论文列表朴素 | 纯文字列表，无缩略图/年份分组/引用数 | 浏览体验差，成果冲击力不足 |
| CV 是一堵"文字墙" | 全部为 markdown 长列表 | 阅读疲劳，重点不突出 |
| 站点描述占位符 | `description: "personal description"` | SEO 与分享预览很差 |
| SEO/社交预览未配置 | `og_image`、`description`、`social.links` 均空 | 微信/Twitter/LinkedIn 分享无缩略图 |
| 大量硬编码 `http://` 绝对链接 | `http://Liuchao-JIN.github.io/files/...` | 非 https、迁移困难、有混合内容风险 |
| "Follow" 按钮无效 | 模板默认残留 | 显得不专业 |
| 无深色模式 / 现代交互 | 默认浅色 | 观感偏旧 |

---

## 二、优先级总览（先做哪些）

**P0（半天内、收益最大）**
1. 修好 `_config.yml` 的 `description` 与 SEO/社交字段（[第九节](#九seo-与可被搜索性)）。
2. 首页 About 增加"研究方向卡片 + 精选成果 + 亮点数字"（[第三节](#三首页about改造)）。
3. 新增 News 动态板块（[第四节](#四新增-news--动态板块)）。
4. 移除/替换无效的 "Follow" 按钮（[第七节](#七视觉与美观主题--配色--排版)）。

**P1（1–2 天）**
5. Publications 加缩略图 + 按年份分组 + Google Scholar 引用徽章（[第五节](#五publications-页面升级)）。
6. CV 页面结构化（时间线 / 分栏 / 折叠长列表）（[第六节](#六cv-页面结构化)）。
7. 配色与字体现代化、深色模式（[第七节](#七视觉与美观主题--配色--排版)）。

**P2（有空再做）**
8. 链接规范化为相对路径 + https（[第十一节](#十一技术债与链接规范)）。
9. 性能与可访问性优化（[第十节](#十性能与可访问性)）。
10. 加分项：访客地图、动画、Publications 可视化等（[第十二节](#十二加分项锦上添花)）。

---

## 三、首页（About）改造

**现状**：`_pages/about.md` 只有一段自我介绍 + 一张封面图，信息密度过低。

**建议结构（从上到下）**：

1. **一句话定位（Hero 标题）**
   - 例如："PhD Candidate @ CUHK · 3D/4D Printing · Smart Materials · Machine-Learning-Driven Design"。
2. **自我介绍**（保留现有段落，稍作润色，突出 HKPFS Fellow 身份）。
3. **研究亮点数字条（Impact numbers）**——用一行醒目的统计增强说服力：

   ```markdown
   <div class="feature__wrapper">
     <div class="feature__item">
       <div class="archive__item">
         <div class="archive__item-body">
           <h2>13+</h2><p>Publications</p>
         </div>
       </div>
     </div>
     <!-- 复制成 4 个：Publications / Citations / Patents / Awards -->
   </div>
   ```
   或更简单地用一行 badge：`📄 13 Papers · 🔬 2 Patents · 🏆 20+ Awards · 📝 30+ Journal Reviews`。

4. **研究方向卡片（Research Interests）**——把当前一句话罗列的方向做成 3–4 张卡片，每张配一句话 + 小图标：
   - 3D / 4D Printing
   - Smart Materials & Metamaterials
   - Soft Robotics
   - ML-Driven Design for Energy Absorption / Harvesting

   可用 Minimal Mistakes 自带的 `feature_row`：
   ```markdown
   {% include feature_row %}
   ```
   配合 About 页 front matter 里定义 `feature_row` 数据。

5. **精选论文（Selected Publications）**——手动挑 3–5 篇代表作，配封面缩略图 + 一句话贡献，放在首页（而不是让访客点进 Publications 才看到）。

6. **最新动态 News**（见[第四节](#四新增-news--动态板块)），放首页折叠展示最近 5 条。

7. **保留封面图**，但建议：
   - 加 `loading="lazy"`；
   - 加图注说明这是哪篇论文的 cover（现在只是 `![Cover Image]`，访客不知含义）；
   - 控制最大宽度，避免在大屏上过大。

**润色示例（About 开头）**：
```markdown
Hi! I'm **Liuchao Jin (金流超)**, a **Hong Kong PhD Fellow** and Ph.D. candidate in
Mechanical & Automation Engineering at **CUHK**, advised by Prof. Wei-Hsin Liao.
I build **machine-learning-driven design methods for 3D/4D printing, smart materials,
metamaterials, and soft robotics** — with applications in energy absorption and harvesting.
```

---

## 四、新增 News / 动态板块

学术主页几乎必备。做法（二选一）：

**方案 A：简单静态列表（推荐先用）**
在 `about.md` 里加一节：
```markdown
## News

- **2025.07** — Named Honorary Alumni Mentor of SCUPI. 🎓
- **2025.05** — Received the PhD IMPAC Award. 🏆
- **2025.02** — Outstanding Students Award at CUHK.
- **2024.08** — New paper in *Applied Materials Today* (ML-driven 4D printing). 📄
- **2024.07** — Started visiting scholar at Shenzhen University.
```
- 用**倒序**（最新在上）；每条尽量带链接（论文/证书/新闻）。
- 建议加"滚动容器 + 固定高度"，超过部分内部滚动，保持首页整洁：
  ```html
  <div style="max-height:260px; overflow-y:auto; padding-right:8px;">
    <!-- news 列表 -->
  </div>
  ```

**方案 B：数据驱动（可维护性更好）**
在 `_data/news.yml` 里维护：
```yaml
- date: "2025-07"
  text: "Named Honorary Alumni Mentor of SCUPI."
  link: "/files/talk_pre/scupi_alumni_mentor.png"
```
再在页面用 Liquid 循环渲染。以后只改数据文件即可。

---

## 五、Publications 页面升级

**现状**：`publications.md` 直接 `for` 循环纯文字，`archive-single.html` 只输出标题+期刊+摘要，无缩略图、无分组。

**建议**：

1. **按年份倒序分组**，加年份小标题：
   ```liquid
   {% assign pubs = site.publications | sort: 'date' | reverse %}
   {% for post in pubs %}
     {% assign year = post.date | date: "%Y" %}
     {% if year != current_year %}
       {% assign current_year = year %}
       <h2 class="pub-year">{{ year }}</h2>
     {% endif %}
     {% include archive-single.html %}
   {% endfor %}
   ```

2. **加论文缩略图（teaser）**——你已经有很多 cover 图（如 `jin2024machine_cover.jpg`）。在每篇 publication 的 front matter 里加：
   ```yaml
   header:
     teaser: /files/my_essay/jin2024machine_cover.jpg
   ```
   `archive-single.html` 已支持 teaser 缩略图渲染，左图右文，观感立刻提升一个档次。

3. **Google Scholar 引用徽章**——用第三方服务（如 [citations.js / Dimensions badge](https://badge.dimensions.ai/)）或手动标注高被引论文：
   `**Cited 30+ times** · [Google Scholar](...)`。

4. **区分类型**：Journal / Conference / Under Review，用小标签（badge）区分，便于快速筛选。

5. **每篇加"一句话贡献"**（TL;DR），比纯摘要更友好：
   ```yaml
   tldr: "First ML framework for 4D-printing inverse design of arbitrary shapes (loss < 0.02 mm)."
   ```

6. **首页只放精选**，完整列表在本页；顶部保留现有"Google Scholar"引导语。

---

## 六、CV 页面结构化

**现状**：`cv.md` 是一整页超长 markdown 列表（教育、经历、论文、专利、奖项 20+、评审 30+、社会服务……），阅读负担重。

**建议**：

1. **顶部加醒目的下载按钮**（现在只是普通链接文字）：
   ```markdown
   <a href="/files/affairs/cv_liuchao_jin.pdf" class="btn btn--primary">📄 Download Full CV (PDF)</a>
   ```

2. **Education / Experience 用时间线样式**（左侧时间轴 + 圆点），比项目符号更直观。可加一小段 CSS `.timeline`。

3. **超长列表折叠**（Awards 20+ 条、Reviewer 30+ 家期刊）——用 `<details>` 折叠，默认只显示前几条：
   ```html
   <details>
     <summary>Show all 30+ journals I review for ▾</summary>
     <!-- 完整列表 -->
   </details>
   ```
   既保留信息量，又不淹没重点。

4. **Honors & Awards 高亮 Top 级奖项**（HKPFS、National Scholarship、Best Paper）用加粗/徽章突出。

5. **加图标分节**：🎓 Education · 💼 Experience · 📄 Publications · 🏆 Awards · 🤝 Service，扫读性更好。

6. 把 `Reviewer` 那一长串期刊改成**多列 / 标签云**展示，而不是竖直长列表。

---

## 七、视觉与美观（主题 / 配色 / 排版）

1. **换配色皮肤**：Minimal Mistakes 内置多套 skin。在 `_config.yml` 加：
   ```yaml
   minimal_mistakes_skin: "default"   # 可选 air / contrast / dirt / neon / plum / sunrise 等
   ```
   建议 `air`（清爽）或自定义一套与 CUHK 紫（#7A003C）呼应的主色，学术又有辨识度。

2. **深色模式**：加一个基于 `prefers-color-scheme` 的 CSS，或引入切换按钮。现代访客好感度高。

3. **字体升级**：正文用更易读的无衬线（如 Inter / Source Han Sans 中英混排），标题可用衬线增强学术感。通过 `_includes/head/custom.html` 引入 Google Fonts。

4. **移除无效 "Follow" 按钮**：`_includes/author-profile.html` 第 23 行
   ```html
   <button class="btn btn--inverse">Follow</button>
   ```
   → 删除，或替换为"📄 Download CV" / "✉ Email me" / "📅 Book with me"（你已有 `bookwithme` 链接）等**真正可用**的行动按钮。

5. **头像**：加圆形边框 + 轻微阴影，`profile_1.png` 建议裁成正方形保证圆形头像不变形。

6. **首页封面图**加圆角与阴影，统一视觉语言。

7. **hover 微交互**：卡片/论文条目加 `transition` + 轻微上浮阴影，提升"精致感"。

---

## 八、导航与信息架构

**现状导航**：CV · Publications · Talks · Teaching · Blog Posts · Amateur。

**建议**：
1. **确认 Talks 是否有内容**——`_config.yml` 里 talks collection 被注释掉了，若无内容建议从导航移除，避免空页。
2. **Portfolio 已有 3 个条目**（PhD/Undergrad 课程资料、文献综述）却被注释隐藏，可考虑改名为 **"Resources"** 或 **"Notes"** 后放出，增加内容丰富度。
3. **加 "News" 锚点**或独立页，方便直达。
4. 导航项建议顺序：**About(Home) · News · Publications · CV · Teaching · Resources · Amateur**，把最重要的 Publications 前置。
5. 移动端确认汉堡菜单正常（模板默认支持）。

---

## 九、SEO 与可被搜索性

**现状**：`description` 是占位符，`og_image`、`social.links`、各 `*_site_verification` 全空 → 分享无预览、搜索排名弱。

**必改（`_config.yml`）**：
```yaml
description: "Liuchao Jin (金流超) — Hong Kong PhD Fellow at CUHK. Research on 3D/4D printing, smart materials, metamaterials, soft robotics, and machine-learning-driven design."

og_image: "/images/profile_1.png"   # 或一张代表性研究图

social:
  type: Person
  name: "Liuchao Jin"
  links:
    - "https://scholar.google.com/citations?user=iYBWir4AAAAJ&hl=en"
    - "https://orcid.org/0000-0001-6204-1922"
    - "https://www.researchgate.net/profile/Liuchao-Jin"
    - "https://www.linkedin.com/in/liuchaojin"
```
**建议补充**：
- 到 [Google Search Console](https://search.google.com/search-console) 验证站点，把值填到 `google_site_verification`。
- 确认 `jekyll-sitemap`、`jekyll-seo-tag` 生效（seo.html 已 include）。
- 每个页面写好 `excerpt`，作为搜索摘要。
- `title` 从 "Jin's Homepage" 改为 **"Liuchao Jin | CUHK"** 之类含全名+机构，利于搜名字。

---

## 十、性能与可访问性

1. **图片**
   - 大图（cover、封面）压缩并转 **WebP**，可显著减小体积。
   - 统一加 `loading="lazy"`。
   - 所有 `<img>` 补充有意义的 `alt`（现在封面 alt 只是 "Cover Image"）。
2. **MathJax**：`head/custom.html` 仍在用 **MathJax 2.7.4**，且无条件全站加载。建议：
   - 升级到 **MathJax 3**（更快）；
   - 仅在需要公式的页面加载。
3. **可访问性（a11y）**
   - 检查配色对比度（WCAG AA）。
   - 链接文字避免"here"，用描述性文字。
   - 键盘可聚焦、focus 样式可见。
4. **Lighthouse 跑一遍**（Chrome DevTools），按报告优化 Performance / SEO / Accessibility 分数。

---

## 十一、技术债与链接规范

**现状**：CV 等页面大量使用硬编码绝对链接：
```
http://Liuchao-JIN.github.io/files/affairs/cv_liuchao_jin.pdf
```
**问题**：① `http://` 非加密（GitHub Pages 强制 https，会有混合内容/跳转）；② 硬编码域名，迁移或本地预览易断。

**建议**：统一改为**站内相对路径**：
```markdown
[Download CV](/files/affairs/cv_liuchao_jin.pdf)
```
可用批量替换（在仓库根目录）：
- 把 `http://Liuchao-JIN.github.io/` 和 `https://liuchao-jin.github.io/` 统一替换为 `/`（注意保留真正的外链）。

**其他**：
- `_config.yml` 里 `analytics.provider: "google-universal"` 但 `tracking_id` 为空 → 要么填 GA4 ID，要么把 provider 设为 `false`，避免加载无效脚本。
- 清理未用到的示例页（`markdown.md`、`terms.md`、`archive-layout-with-content.md` 等模板残留），保持仓库整洁。

---

## 十二、加分项（锦上添花）

- **访客足迹地图 / 访问统计**：如 [ClustrMaps](https://clustrmaps.com/) 或 [Mapmyvisitors](https://mapmyvisitors.com/)，放侧栏，直观展示影响力。
- **Publications 可视化**：一张按年份的论文数量柱状图，或研究关键词词云。
- **"News" 时间线动画**：滚动淡入。
- **Selected Media / 视频嵌入**：Amateur 页把 YouTube/Bilibili 视频**内嵌播放**（而非纯链接），观感更好。
- **多语言**：加中文版首页（你有中英双重受众：CUHK + SCU/SUSTech/SZU）。
- **"Now" 页面**：一句话说明"我目前在做什么"，增加人情味。
- **RSS/订阅**：Blog 已装 `jekyll-feed`，在页脚放订阅入口。
- **Favicon/主题色**：`theme-color` 已是白色，可改为品牌主色，移动端地址栏更统一。

---

## 十三、落地清单（Checklist）

复制到你的 issue / todo 里逐项打勾：

**P0**
- [ ] `_config.yml`：填写真实 `description`
- [ ] `_config.yml`：配置 `og_image` + `social.links`
- [ ] `_config.yml`：`title` 改为含全名+机构
- [ ] `_config.yml`：`analytics.tracking_id` 填写或将 provider 设为 false
- [ ] About 首页：加研究方向卡片 + 亮点数字 + 精选论文
- [ ] 新增 News 板块（首页或独立页）
- [ ] 删除无效 "Follow" 按钮，替换为 CV / Email / Book 按钮

**P1**
- [ ] Publications：加 teaser 缩略图（各 publication front matter）
- [ ] Publications：按年份分组 + 类型标签
- [ ] CV：加下载按钮 + 时间线 + 长列表折叠(`<details>`)
- [ ] 选定 `minimal_mistakes_skin` 皮肤 / 自定义配色
- [ ] 字体升级 + 深色模式

**P2**
- [ ] 链接统一改为相对路径 + https
- [ ] 图片压缩 / WebP / lazy-load / alt
- [ ] MathJax 升级到 v3 并按需加载
- [ ] 跑 Lighthouse 并按报告优化
- [ ] 清理模板残留页面
- [ ] 加分项：访客地图 / 视频内嵌 / 中文版 等

---

### 备注

- 本站基于 Jekyll，**本地预览**：`bundle exec jekyll serve`，浏览器打开 `http://localhost:4000`，改动实时可见（改 `_config.yml` 需重启）。
- 所有改动建议**先在分支上做**，确认无误再合并到 `master`/`main`，避免线上主页短暂异常。
- 如果需要，我可以按这份方案**逐项帮你落地实现**（例如先从 P0 的首页改造 + News + SEO 三项开始）。
